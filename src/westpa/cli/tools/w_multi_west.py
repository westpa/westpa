import logging
import numpy as np
import pickle

log = logging.getLogger(__name__)
from westpa.tools.core import WESTTool
from westpa.core.data_manager import n_iter_dtype, istate_dtype
from westpa.tools.progress import ProgressIndicatorComponent
from westpa.core import h5io
from westpa.tools.core import WESTMultiTool

# from westtools.dtypes import iter_block_ci_dtype as ci_dtype
import gc

# from pympler.tracker import SummaryTracker

ci_dtype = np.dtype(
    [
        ('iter_start', n_iter_dtype),
        ('iter_stop', n_iter_dtype),
        ('expected', np.float64),
        ('ci_lbound', np.float64),
        ('ci_ubound', np.float64),
        ('corr_len', n_iter_dtype),
        ('variance', np.float64),
        ('stderrormean', np.float64),
    ]
)

# directory locations are stored in a .yaml file with this format:
# ---
# PATHS: ['/path/to/simulation/1','/path/to/simulation/2',...,
# '/path/to/simulation/n']

# Straight up stolen from the data manager.  In the future, maybe I can just sort it by subbing in the appropriate values.


def get_bin_mapper(we_h5file, hashval):
    '''Look up the given hash value in the binning table, unpickling and returning the corresponding
    bin mapper if available, or raising KeyError if not.'''

    # Convert to a hex digest if we need to
    try:
        hashval = hashval.hexdigest()
    except AttributeError:
        pass

    while True:
        # these will raise KeyError if the group doesn't exist, which also means
        # that bin data is not available, so no special treatment here
        try:
            binning_group = we_h5file['/bin_topologies']
            index = binning_group['index']
            pkl = binning_group['pickles']
        except KeyError:
            raise KeyError('hash {} not found. Could not retrieve binning group'.format(hashval))

        n_entries = len(index)
        if n_entries == 0:
            raise KeyError('hash {} not found. No entries in index'.format(hashval))

        chunksize = 1024

        for istart in range(0, n_entries, chunksize):
            chunk = index[istart : min(istart + chunksize, n_entries)]
            for i in range(len(chunk)):
                if chunk[i]['hash'] == hashval:
                    pkldat = bytes(pkl[istart + i, 0 : chunk[i]['pickle_len']].data)
                    mapper = pickle.loads(pkldat)
                    log.debug('loaded {!r} from {!r}'.format(mapper, binning_group))
                    log.debug('hash value {!r}'.format(hashval))
                    return mapper

        raise KeyError('hash {} not found'.format(hashval))


def create_idtype_array(input_array):
    '''Return a new array with the new istate_dtype while preserving old data.'''
    new_array = np.zeros(input_array.shape, dtype=istate_dtype)
    for j in input_array.dtype.names:
        new_array[j] = input_array[j].copy()

    # Need to turn 'basis_auxref' to empty bytestrings...
    new_array['basis_auxref'] = b''

    return new_array


class WMultiWest(WESTMultiTool):
    prog = 'w_multi_west'
    description = '''\
Tool designed to combine multiple WESTPA simulations while accounting for
reweighting.
-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WESTTool, self).__init__()
        self.progress = ProgressIndicatorComponent()
        # We no longer care about a lot of this.
        self.ntrials = 0
        self.nstates = 0
        self.kin_trial = {}
        self.west = {}
        self.niters = 0

    def add_args(self, parser):
        self.progress.add_args(parser)
        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('-o', '--output-file', default='multi.h5', help='''The name of the output file to store results in.''')
        iogroup.add_argument(
            '-W',
            '--west',
            '--WEST_H5FILE',
            default='west.h5',
            help='''The name of the main .h5 file inside each simulation
                             directory''',
        )
        iogroup.add_argument('-a', '--aux', action='append', help='''Names of additional auxiliary datasets to be combined''')
        iogroup.add_argument('-aa', '--auxall', action='store_true', help='''Combine all auxiliary datasets. Default: False''')
        iogroup.add_argument('-nr', '--no-reweight', action='store_true', help='''Do not reweight. Default: False''')
        iogroup.add_argument(
            '-ib', '--ibstates', action='store_true', help='''Attempt to combine ibstates dataset. Default: False'''
        )

    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_file, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)

        opened_files = self.generate_file_list([self.west])
        self.westH5 = opened_files[self.west]
        # Just some temp things while I clean everything up...
        # west_files = self.westH5
        # Determine max iteration ...

        # We can't really use the old method anymore, as we need to calculate rates in the bootstrap.
        # Ergo, we're going to load things like w_kinavg, but that's all.
        # We'll just load them up and store them internally, for the moment.

    def process_args(self, args):
        self.progress.process_args(args)
        self.output_file = args.output_file
        self.output_file_name = args.output_file
        self.west = args.west
        self.sims = args.sims
        self.aux = args.aux
        self.auxall = args.auxall
        self.reweight = args.no_reweight
        self.ibstates = args.ibstates

    def total_number_of_walkers(self):
        self.total_walkers = [0] * self.niters
        for key, west in self.westH5.items():
            # Sometimes, we're smaller or larger by one.  Hm.
            try:
                self.total_walkers[:] += west['summary'][:-1]['n_particles']
            except (ValueError):
                self.total_walkers[:] += west['summary'][:-1]['n_particles'][: len(self.total_walkers)]

    class Segment:
        def __init__(self, weight=0, iteration=0, simid=0, recycled_in=0):
            self.weight = weight
            self.iteration = iteration
            self.simid = simid
            self.recycled_in = recycled_in

    def go(self):
        pi = self.progress.indicator
        self.istates = True  # Assume serendipitously istates is same between runs...
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            self.total_number_of_walkers()
            if self.auxall is True:
                self.aux = list(self.westH5[1]['iterations/iter_00000001/auxdata'].keys())
            # Create a giant WEST.h5 file, separating the individual walkers, and renormalizing the weights.
            # It should then be compatible with existing toolsets.
            # Isn't really going to start with auxdata, but we'll add it in.

            # self.niters = 500
            # Initialize data manager...
            # Just bullshit for the current system.
            # self.niters = self.westH5[1].attrs['west_current_iteration'] - 1
            # print(self.niters, len(self.westH5))
            # self.data_manager = data_manager.WESTDataManager()
            westh5 = []
            self.source_sinks = []
            self.n_sims = {}
            istate_addition = [0]
            for ifile, (key, west) in enumerate(self.westH5.items()):
                d = {'west': west, 'wm': None, 'rt': None, 'remove_next_cycle': [], 'seg_index': None}
                # We're getting the bin mapper, then setting the recycling target...
                binhash = west['iterations/iter_{0:08d}'.format(2)].attrs['binhash']
                bin_mapper = get_bin_mapper(west, bytes(binhash, 'utf-8'))
                try:
                    d['rt'] = bin_mapper.assign(west['tstates']['0']['pcoord'][...])[0]
                    self.source_sinks.append(bin_mapper.assign(west['tstates']['0']['pcoord'][...])[0])
                except KeyError:
                    d['rt'] = None
                    self.source_sinks.append(None)
                    pass
                # We're going to make a set of source and sink states that we can iterate through, eventually.
                # Keep a count of how many simulations for this particular recycling target we have...
                try:
                    self.n_sims[d['rt']] += 1
                except KeyError:
                    self.n_sims[d['rt']] = 1
                westh5.append(d)
                if ifile == 0:
                    self.niters = west.attrs['west_current_iteration'] - 1
                else:
                    self.niters = min(west.attrs['west_current_iteration'] - 1, self.niters)

                istate_addition.append(istate_addition[-1] + len(west['ibstates/0/istate_index']))
                # Check to see if all the bstates are identical
                if self.ibstates:
                    check = [False, False]  # Assuming they're false, so not accidentally outputing anything that errors out.
                    try:
                        check[0] = np.array_equal(bstate_index, west['ibstates/0/bstate_index'][:])
                        check[1] = np.array_equal(bstate_pcoord, west['ibstates/0/bstate_pcoord'][:])
                        if not np.all(check):
                            print(f'File {ifile} used different bstates than the first file. Will skip exporting ibstates dataset.')
                            self.ibstates = False
                    except NameError:
                        bstate_index = west['ibstates/0/bstate_index'][:]  # noqa: F841
                        bstate_pcoord = west['ibstates/0/bstate_pcoord'][:]  # noqa: F841

            start_point = []
            self.source_sinks = list(set(self.source_sinks))
            # We'll need a global list of walkers to add to and take care of during the next round of simulations, as well as the current one.
            # We'll organize it by source and sink states.
            self.past_iter = {}
            self.futr_iter = {}
            self.past_rm = {}
            self.futr_rm = {}
            for i in self.source_sinks:
                self.past_iter[i] = []
                self.futr_iter[i] = []
                self.past_rm[i] = []
                self.futr_rm[i] = []
            print(pi.new_operation('Recreating...', self.niters))
            # tracker = SummaryTracker()
            # self.output_file.close()

            if self.ibstates:
                # Copying the ibstates group from the first file as base
                self.output_file.copy(self.westH5[1]['ibstates'], self.output_file)
                del self.output_file['ibstates/0/istate_pcoord']
                del self.output_file['ibstates/0/istate_index']

                # Combining the rest of the istate datasets
                for ifile, (key, west) in enumerate(self.westH5.items()):
                    if ifile == 0:
                        final_istate_index = west['ibstates/0/istate_index']
                        final_istate_pcoord = west['ibstates/0/istate_pcoord']
                        if final_istate_index.dtype != istate_dtype:
                            final_istate_index = create_idtype_array(final_istate_index)
                    else:
                        addition = west['ibstates/0/istate_index'][:]
                        if addition.dtype != istate_dtype:
                            addition = create_idtype_array(addition)
                        final_istate_index = np.append(final_istate_index, addition)
                        final_istate_pcoord = np.append(final_istate_pcoord, west['ibstates/0/istate_pcoord'][:])

                # Saving them into self.output_file
                self.output_file['ibstates/0'].create_dataset('istate_index', data=final_istate_index, dtype=istate_dtype)
                self.output_file['ibstates/0'].create_dataset('istate_pcoord', data=final_istate_pcoord)

            for iter in range(self.niters):
                # We have the following datasets in each iteration:
                # ibstates, which can now be combined with --ibstates
                # pcoord
                # seg_index
                # wtgraph
                # wtgraph is going to be a little more complex to handle, but not too bad.
                # aux data specified
                iter += 1
                ifile = 0
                # self.output_file = h5io.WESTPAH5File(self.output_file_name, 'w', creating_program=True)
                # Determine how many simulations to append or remove per west file.
                # self.segments = {}
                # for key,value in self.n_sims.items():
                #    self.segments[key] = int(np.floor(len(self.past_iter[key]) / value))

                # run_once = 0
                # total_current_sims = 0
                # for i in self.source_sinks:
                #    total_current_sims += len(self.past_iter[i])
                #    total_current_sims += len(self.past_rm[i])
                for ifile, west in enumerate(westh5):
                    westdict = west['west']
                    seg_index = westdict['iterations/iter_{0:08d}'.format(iter)]['seg_index'][...]
                    pcoord = westdict['iterations/iter_{0:08d}'.format(iter)]['pcoord'][...]
                    wtgraph = westdict['iterations/iter_{0:08d}'.format(iter)]['wtgraph'][...]
                    # new_weight = westdict['iterations/iter_{0:08d}'.format(iter)]['new_weight'][...]
                    if self.aux:
                        auxdata = {}
                        for i in self.aux:
                            auxdata[str(i)] = westdict['iterations/iter_{0:08d}'.format(iter)]['auxdata'][str(i)][...]
                    if iter == 1 and ifile == 0:
                        new_dtype = np.dtype(seg_index.dtype.descr + [('group', '<i8')])
                    new_seg_index = np.zeros(seg_index.shape, dtype=new_dtype)
                    for dt, val in seg_index.dtype.fields.items():
                        new_seg_index[dt] = seg_index[dt]
                    new_seg_index['group'] = ifile
                    del seg_index
                    seg_index = new_seg_index[...]
                    del new_seg_index
                    if ifile == 0:
                        mseg = seg_index
                        mpco = pcoord
                        mwtg = wtgraph
                        if self.aux:
                            maux = {}
                            for i in self.aux:
                                maux[str(i)] = auxdata[str(i)]
                        if iter == 1:
                            summary = westdict['summary'][...]

                        start_point.append(0)
                    if ifile != 0:
                        # print(mseg.shape, seg_index.shape, ifile)
                        # print(mpco.shape, pcoord.shape, ifile)
                        # print(mwtg.shape, wtgraph.shape, ifile)
                        if iter != 1:
                            addition = prev_start_point[ifile]  # noqa: F821
                        else:
                            addition = mseg.shape[0]
                        seg_index['parent_id'][np.where(seg_index['parent_id'] >= 0)] += addition
                        seg_index['parent_id'][np.where(seg_index['parent_id'] < 0)] -= istate_addition[ifile]
                        seg_index['wtg_offset'] += mwtg.shape[0]
                        start_point.append(mseg.shape[0])
                        wtgraph += mwtg.shape[0]
                        mseg = np.concatenate((mseg, seg_index))
                        mpco = np.concatenate((mpco, pcoord))
                        mwtg = np.concatenate((mwtg, wtgraph))
                        if self.aux:
                            for i in self.aux:
                                maux[str(i)] = np.concatenate((maux[str(i)], auxdata[str(i)]))
                    ifile += 1
                    del seg_index, pcoord, wtgraph, westdict
                    if self.aux:
                        del auxdata
                gc.collect()
                # Make a real copy to use in the next iteration.
                # self.past_iter = self.futr_iter.copy()
                # self.past_rm[i] = self.futr_rm.copy()
                prev_start_point = start_point  # noqa: F841
                start_point = []
                # This is... maybe wrong, actually?  Or at least, it's not ALL that is required for normalizing things.
                # We need to weight everything by 1/N, then just normalize if that normalization was wrong.  Keep the relative weights sane.
                # ... or actually, no, that's fine, nevermind, what's wrong with me?  But we'll leave it in for now.

                # Normalize weight of each iteration, done unless specified not to.
                if not self.reweight:
                    mseg['weight'] /= mseg['weight'].sum()

                summary['n_particles'][iter - 1] = mseg.shape[0]
                summary['norm'][iter - 1] = mseg['weight'].sum()
                summary['min_seg_prob'][iter - 1] = min(mseg['weight'])
                summary['max_seg_prob'][iter - 1] = max(mseg['weight'])

                curr_iter = self.output_file.create_group('iterations/iter_{0:08d}'.format(iter))
                curr_iter.attrs['n_iter'] = iter

                # Hard-link ibstates dataset to the main one
                if self.ibstates:
                    curr_iter['ibstates'] = self.output_file['ibstates/0']

                ds_rate_evol = curr_iter.create_dataset('wtgraph', data=mwtg, shuffle=True, compression=9)
                ds_rate_evol = curr_iter.create_dataset('seg_index', data=mseg, shuffle=True, compression=9)
                ds_rate_evol = curr_iter.create_dataset('pcoord', data=mpco, shuffle=True, compression=9)
                if self.aux:
                    aux_iter = self.output_file.create_group('iterations/iter_{0:08d}/auxdata'.format(iter))
                    for i in self.aux:
                        ds_rate_evol = aux_iter.create_dataset(str(i), data=maux[str(i)], shuffle=True, compression=9)
                # We need to be careful about memory, here.  We are blowing uppppp.
                # We're STILL blowing up.  Criiiiiipes.
                # self.segments = {}
                del mseg, mpco, mwtg, ds_rate_evol, curr_iter  # , self.segments
                if self.aux:
                    del maux, aux_iter
                gc.collect()
                self.output_file.flush()
                # self.output_file.close()
                # print("How big is our summary?")
                # print(sys.getsizeof(summary))
                # objgraph.show_most_common_types(limit=50)
                # objgraph.show_growth(limit=10)
                # objgraph.show_most_common_types(objects=objgraph.get_leaking_objects())
                pi.progress += 1

        pi.new_operation('Writing to file...')
        ds_rate_evol = self.output_file.create_dataset('summary', data=summary, shuffle=True, compression=9)  # noqa: F841
        self.output_file.attrs['west_current_iteration'] = self.niters
        self.output_file.attrs['west_file_format_version'] = 7
        self.output_file.attrs['west_iter_prec'] = 8
        self.output_file.attrs['westpa_fileformat_version'] = 7
        self.output_file.attrs['westpa_iter_prec'] = 8


def entry_point():
    WMultiWest().main()


if __name__ == '__main__':
    entry_point()
