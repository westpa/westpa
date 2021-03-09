import h5py
import numpy as np

from westpa.tools import WESTTool, WESTDataReader, IterRangeSelection, ProgressIndicatorComponent

import westpa
from westpa.core import h5io
from westpa.core.data_manager import seg_id_dtype, n_iter_dtype, weight_dtype
from westpa.core.binning import assignments_list_to_table


class WNTopTool(WESTTool):
    prog = 'w_ntop'
    description = '''\
Select walkers from bins . An assignment file mapping walkers to
bins at each timepoint is required (see``w_assign --help`` for further
information on generating this file). By default, high-weight walkers are
selected (hence the name ``w_ntop``: select the N top-weighted walkers from
each bin); however, minimum weight walkers and randomly-selected walkers
may be selected instead.


-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "ntop.h5") contains the following
datasets:

  ``/n_iter`` [iteration]
    *(Integer)* Iteration numbers for each entry in other datasets.

  ``/n_segs`` [iteration][bin]
    *(Integer)* Number of segments in each bin/state in the given iteration.
    This will generally be the same as the number requested with
    ``--n/--count`` but may be smaller if the requested number of walkers
    does not exist.

  ``/seg_ids`` [iteration][bin][segment]
    *(Integer)* Matching segments in each iteration for each bin.
    For an iteration ``n_iter``, only the first ``n_iter`` entries are
    valid. For example, the full list of matching seg_ids in bin 0 in the
    first stored iteration is ``seg_ids[0][0][:n_segs[0]]``.

  ``/weights`` [iteration][bin][segment]
    *(Floating-point)* Weights for each matching segment in ``/seg_ids``.


-----------------------------------------------------------------------------
Command-line arguments
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super().__init__()

        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.progress = ProgressIndicatorComponent()
        self.output_file = None
        self.assignments_filename = None
        self.output_filename = None
        self.what = None
        self.timepoint = None
        self.count = None

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)

        igroup = parser.add_argument_group('input options')
        igroup.add_argument(
            '-a',
            '--assignments',
            default='assign.h5',
            help='''Use assignments from the given ASSIGNMENTS file (default: %(default)s).''',
        )

        sgroup = parser.add_argument_group('selection options')
        sgroup.add_argument(
            '-n',
            '--count',
            type=int,
            default=1,
            help='''Select COUNT walkers from each iteration for each bin (default: %(default)s).''',
        )
        sgroup.add_argument(
            '-t',
            '--timepoint',
            type=int,
            default=-1,
            help='''Base selection on the given TIMEPOINT within each iteration. Default (-1)
                            corresponds to the last timepoint.''',
        )
        cgroup = parser.add_mutually_exclusive_group()
        cgroup.add_argument(
            '--highweight',
            dest='select_what',
            action='store_const',
            const='highweight',
            help='''Select COUNT highest-weight walkers from each bin.''',
        )
        cgroup.add_argument(
            '--lowweight',
            dest='select_what',
            action='store_const',
            const='lowweight',
            help='''Select COUNT lowest-weight walkers from each bin.''',
        )
        cgroup.add_argument(
            '--random',
            dest='select_what',
            action='store_const',
            const='random',
            help='''Select COUNT walkers randomly from each bin.''',
        )
        parser.set_defaults(select_what='highweight')

        ogroup = parser.add_argument_group('output options')
        ogroup.add_argument('-o', '--output', default='ntop.h5', help='''Write output to OUTPUT (default: %(default)s).''')
        self.progress.add_args(parser)

    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)
        self.what = args.select_what
        self.output_filename = args.output
        self.assignments_filename = args.assignments
        self.count = args.count
        self.timepoint = args.timepoint

    def go(self):
        self.data_reader.open('r')
        assignments_file = h5py.File(self.assignments_filename, mode='r')
        output_file = h5io.WESTPAH5File(self.output_filename, mode='w')
        pi = self.progress.indicator
        count = self.count
        timepoint = self.timepoint

        nbins = assignments_file.attrs['nbins'] + 1
        assignments_ds = assignments_file['assignments']

        iter_start, iter_stop = self.iter_range.iter_start, self.iter_range.iter_stop
        iter_count = iter_stop - iter_start
        h5io.check_iter_range_least(assignments_ds, iter_start, iter_stop)
        nsegs = assignments_file['nsegs'][h5io.get_iteration_slice(assignments_file['nsegs'], iter_start, iter_stop)]

        output_file.create_dataset('n_iter', dtype=n_iter_dtype, data=list(range(iter_start, iter_stop)))

        seg_count_ds = output_file.create_dataset('nsegs', dtype=np.uint, shape=(iter_count, nbins))
        matching_segs_ds = output_file.create_dataset(
            'seg_ids',
            shape=(iter_count, nbins, count),
            dtype=seg_id_dtype,
            chunks=h5io.calc_chunksize((iter_count, nbins, count), seg_id_dtype),
            shuffle=True,
            compression=9,
        )
        weights_ds = output_file.create_dataset(
            'weights',
            shape=(iter_count, nbins, count),
            dtype=weight_dtype,
            chunks=h5io.calc_chunksize((iter_count, nbins, count), weight_dtype),
            shuffle=True,
            compression=9,
        )
        what = self.what

        with pi:
            pi.new_operation('Finding matching segments', extent=iter_count)
            for iiter, n_iter in enumerate(range(iter_start, iter_stop)):
                assignments = np.require(
                    assignments_ds[h5io.get_iteration_entry(assignments_ds, n_iter) + np.index_exp[:, timepoint]],
                    dtype=westpa.core.binning.index_dtype,
                )
                all_weights = self.data_reader.get_iter_group(n_iter)['seg_index']['weight']

                # the following Cython function just executes this loop:
                # for iseg in xrange(nsegs[iiter]):
                #    segs_by_bin[iseg,assignments[iseg]] = True
                segs_by_bin = assignments_list_to_table(nsegs[iiter], nbins, assignments)
                for ibin in range(nbins):
                    segs = np.nonzero(segs_by_bin[:, ibin])[0]

                    seg_count_ds[iiter, ibin] = min(len(segs), count)

                    if len(segs):
                        weights = all_weights.take(segs)

                        if what == 'lowweight':
                            indices = np.argsort(weights)[:count]
                        elif what == 'highweight':
                            indices = np.argsort(weights)[::-1][:count]
                        else:
                            assert what == 'random'
                            indices = np.random.permutation(len(weights))

                        matching_segs_ds[iiter, ibin, : len(segs)] = segs.take(indices)
                        weights_ds[iiter, ibin, : len(segs)] = weights.take(indices)
                        del segs, weights

                del assignments, segs_by_bin, all_weights
                pi.progress += 1


def entry_point():
    WNTopTool().main()


if __name__ == '__main__':
    entry_point()
