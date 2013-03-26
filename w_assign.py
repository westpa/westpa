from __future__ import print_function, division; __metaclass__ = type
import sys
import logging
import math

from west.data_manager import seg_id_dtype
from westpa.binning.assign import index_dtype, UNKNOWN_INDEX
from westpa.binning._assign import assign_and_label #@UnresolvedImport
from westtools.tool_classes import WESTTool, WESTParallelTool, WESTDataReader, BinMappingComponent
import numpy, h5py
from westpa import h5io
from westpa.h5io import WESTPAH5File
from westpa.extloader import get_object
import westpa

log = logging.getLogger('westtools.w_assign')

def parse_pcoord_value(pc_str):
    namespace = {'math': math,
                 'numpy': numpy,
                 'inf': float('inf')}
    
    arr = numpy.array(eval(pc_str,namespace))
    if arr.ndim == 0:
        arr.shape = (1,1)
    elif arr.ndim == 1:
        arr.shape = (1,) + arr.shape 
    else:
        raise ValueError('too many dimensions')
    return arr

def default_construct_pcoord(n_iter, iter_group):
    return iter_group['pcoord'][...]


def _assign_and_label(nsegs_lb, nsegs_ub, npts, parent_ids,
                      bm, state_map, last_labels, pcoords):
    

    
    assignments, trajlabels = assign_and_label(nsegs_lb, nsegs_ub, npts, parent_ids,
                                               bm.assign, state_map, last_labels, pcoords)

    return (assignments, trajlabels, nsegs_lb, nsegs_ub)

class WAssign(WESTParallelTool):
    prog='w_assign'
    description = '''\
Assign walkers to bins, producing a file (by default named "assign.h5")
which can be used in subsequent analysis.

Progress coordinate data is taken by default from the "pcoord" dataset
for each iteration in the main HDF5 file (usually west.h5). However,
an arbitrary function can be provided with -p/--construct-pcoord,
which can be used to consolidate data from several sources or otherwise
preprocess it before binning occurs. 

For consistency in subsequent analysis operations, the entire dataset
must be assigned, even if only a subset of the data will be used. This
ensures that analyses that rely on tracing trajectories always know the
originating bin of each trajectory.

Optionally, kinetic macrostates may be defined in terms of sets of bins.
Each trajectory will be labeled with the kinetic macrostate it was most
recently in at each timepoint, for use in subsequent kinetic analysis.
Macrostates corresponding to single bins may be identified on the command
line. States corresponding to multiple bins use a YAML input file consisting
of a list of states, each with a name and a list of coordinate tuples; bins
containing these coordinates will be mapped to the containing state. For
instance, the following file::

    ---
    states:
      - label: unbound
        coords:
          - [9.0, 1.0]
          - [9.0, 2.0]
      - label: bound
        coords:
          - [0.1, 0.0]

produces two macrostates: the first state is called "unbound" and consists of
bins containing the (2-dimensional) progress coordinate values (9.0, 1.0) and
(9.0, 2.0); the second state is called "bound" and consists of the single bin
containing the point (0.1, 0.0).
'''
    
    def __init__(self):
        super(WAssign,self).__init__()
        self.data_reader = WESTDataReader() 
        self.binning = BinMappingComponent()
        self.output_file = None
        self.construct_pcoord = default_construct_pcoord
        self.states = []
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        
        self.binning.add_args(parser, suppress=['--bins-from-file'])
        
        agroup = parser.add_argument_group('other options')
        
        agroup.add_argument('-p', '--construct-pcoord', metavar='MODULE.FUNCTION',
                             help='''Use the given function to construct progress coordinate data
                             for each iteration. This function will be called once per iteration as
                             ``get_pcoord(n_iter, iter_group)``, and must return an array indexable as
                             [seg_id][timepoint][dimension]. By default, the "pcoord" dataset is
                             loaded into memory and returned.''')
        
        agroup.add_argument('-o', '--output', dest='output', default='assign.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        
        sgroup = agroup.add_mutually_exclusive_group()
        sgroup.add_argument('--states', nargs='+', metavar='STATEDEF',
                            help='''Single-bin kinetic macrostate, specified by a coordinate tuple (e.g. '1.0' or '[1.0,1.0]'),
                            optionally labeled (e.g. 'bound:[1.0,1.0]'). States corresponding to multiple bins
                            must be specified with --states-from-file.''')
        sgroup.add_argument('--states-from-file', metavar='STATEFILE',
                            help='''Load kinetic macrostates from the YAML file STATEFILE. See description
                            above for the appropriate structure.''')        
        sgroup.add_argument('--states-from-function', metavar='STATEFUNC',
                            help='''Load kinetic macrostates from the function STATEFUNC, specified as
                            module_name.func_name. This function is called with the bin mapper as an argument,
                            and must return a list of dictionaries {'label': state_label, 'coords': 2d_array_like}
                            one for each macrostate; the 'coords' entry must contain enough rows to identify all bins
                            in the macrostate.''')        


    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open(mode='r')
        self.binning.process_args(args)
        if args.construct_pcoord:
            self.construct_pcoord = get_object(args.construct_pcoord,path=['.'])
            
        if args.states:
            self.parse_cmdline_states(args.states)
        elif args.states_from_file:
            self.load_state_file(args.states_from_file)
        elif args.states_from_function:
            self.load_states_from_function(get_object(args.states_from_function,path=['.']))
            
        if self.states and len(self.states) < 2:
            raise ValueError('zero, two, or more macrostates are required')

        self.output_file = WESTPAH5File(args.output, 'w', creating_program=True)
        log.debug('state list: {!r}'.format(self.states))
        
    def parse_cmdline_states(self, state_strings):
        states = []
        for istring, state_string in enumerate(state_strings):
            try:
                (label, coord_str) = state_string.split(':')
            except ValueError:
                label = 'state{}'.format(istring)
                coord_str = state_string
            coord = parse_pcoord_value(coord_str)
            states.append({'label': label, 'coords': coord})
        self.states = states
    
    def load_state_file(self, state_filename):
        import yaml
        ydict = yaml.load(open(state_filename, 'rt'))
        ystates = ydict['states']
        
        states = []
        for istate, ystate in enumerate(ystates):
            state = {}
            state['label'] = ystate.get('label', 'state{}'.format(istate))
            # coords can be:
            #  - a scalar, in which case it is one bin, 1-D
            #  - a single list, which is rejected as ambiguous
            #  - a list of lists, which is a list of coordinate tuples
            coords = numpy.array(ystate['coords'])
            if coords.ndim == 0:
                coords.shape = (1,1)
            elif coords.ndim == 1:
                raise ValueError('list {!r} is ambiguous (list of 1-d coordinates, or single multi-d coordinate?)'
                                 .format(ystate['coords']))
            elif coords.ndim > 2:
                raise ValueError('coordinates must be 2-D')
            state['coords'] = coords
            states.append(state)
        self.states = states

    def assign_iteration(self, nsegs, npts, parent_ids, state_map, last_labels, pcoords, n_workers, n_iter):
        ''' Method to encapsulate the segment slicing (into n_worker slices) and parallel job submission
            Submits job(s), waits on completion, splices them back together
            Returns: assignments, trajlabels for this iteration'''

        futures = []

        #Submit jobs to work manager
        for islice in xrange(n_workers):
            lb = (islice)*(nsegs//n_workers)
            ub = (islice+1)*(nsegs//n_workers)
            futures.append(self.work_manager.submit(_assign_and_label, 
                           args=(lb, ub, npts, parent_ids, self.binning.mapper, state_map, last_labels, pcoords)))
       
        assignments = numpy.empty((nsegs, npts), dtype=index_dtype)
        trajlabels = numpy.empty((nsegs, npts), dtype=index_dtype)

        nrecvd = 0
        for future in self.work_manager.as_completed(futures):
            nrecvd += 1
            assign_slice, traj_slice, lb, ub = future.get_result()

            assignments[lb:ub, :] = assign_slice
            trajlabels[lb:ub, :] = traj_slice

        assert nrecvd == n_workers, "Error: iteration {} did not complete successfully".format(n_iter)

        return (assignments, trajlabels)

        
    def load_states_from_function(self, statefunc):
        states = statefunc(self.binning.mapper)
        for istate, state in enumerate(states):
            state.setdefault('label','state{}'.format(istate))
            try:
                state['coords'] = numpy.array(state['coords'])
            except KeyError:
                raise ValueError('state function {!r} returned a state {!r} without coordinates'.format(statefunc,state))
        self.states = states
        log.debug('loaded states: {!r}'.format(self.states))
        
    def go(self):
        assign = self.binning.mapper.assign
        
        
        iter_start = 1 
        iter_stop =  self.data_reader.current_iteration
        
        h5io.stamp_iter_range(self.output_file, iter_start, iter_stop)
        
        nbins = self.binning.mapper.nbins
        self.output_file.attrs['nbins'] = nbins 
        
        state_map = None
        if self.states:
            state_map = numpy.empty((self.binning.mapper.nbins,), index_dtype)
            state_map[:] = UNKNOWN_INDEX
            
            state_labels = [state['label'] for state in self.states]
            
            for istate, sdict in enumerate(self.states):
                assert state_labels[istate] == sdict['label'] #sanity check
                state_assignments = assign(sdict['coords'])
                for assignment in state_assignments:
                    state_map[assignment] = istate
            
            # Don't compress this if it's tiny
            # It's a microoptimization, but a simple one.
            if nbins <= 10:
                dsopts = {}
            else:
                dsopts = {'compression': 9,
                          'shuffle': True}
            self.output_file.create_dataset('state_map', data=state_map, **dsopts)
            self.output_file['state_labels'] = state_labels
        
        iter_count = iter_stop - iter_start
        nsegs = numpy.empty((iter_count,), seg_id_dtype)
        npts = numpy.empty((iter_count,), seg_id_dtype)
        
        # scan for largest number of segments and largest number of points
        for iiter, n_iter in enumerate(xrange(iter_start,iter_stop)):
            iter_group = self.data_reader.get_iter_group(n_iter)
            nsegs[iiter], npts[iiter] = iter_group['pcoord'].shape[0:2]
            del iter_group
            
        # create datasets
        self.output_file.create_dataset('nsegs', data=nsegs, shuffle=True, compression=9)
        self.output_file.create_dataset('npts', data=npts, shuffle=True, compression=9)
        
        max_nsegs = nsegs.max()
        max_npts = npts.max()
        
        assignments_shape = (iter_count,max_nsegs,max_npts)
        assignments_ds = self.output_file.create_dataset('assignments', dtype=index_dtype, shape=assignments_shape,
                                                         compression=9, shuffle=True,
                                                         chunks=h5io.calc_chunksize(assignments_shape, index_dtype),
                                                         fillvalue=UNKNOWN_INDEX)
        if self.states:
            trajlabels_ds = self.output_file.create_dataset('trajlabels', dtype=index_dtype, shape=assignments_shape,
                                                            compression=9, shuffle=True,
                                                            chunks=h5io.calc_chunksize(assignments_shape, index_dtype),
                                                            fillvalue=UNKNOWN_INDEX)

        last_labels = None # mapping of seg_id to last macrostate inhabited      
        for iiter, n_iter in enumerate(xrange(iter_start,iter_stop)):

            #get iteration info in this block
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rIteration {}'.format(n_iter),end='')
                sys.stdout.flush()
            iter_group = self.data_reader.get_iter_group(n_iter)
            parent_ids = iter_group['seg_index']['parent_id']
            pcoords = self.construct_pcoord(n_iter,iter_group)
            
            assert pcoords.shape[0] == nsegs[iiter]
            assert pcoords.shape[1] == npts[iiter]
            
            if iiter == 0:
                last_labels = numpy.empty((nsegs[iiter],), index_dtype)
                last_labels[:] = UNKNOWN_INDEX

            #Determine the number of available workers
            assert self.work_manager, 'No work manager created for {!r}'.format(self)
            n_workers = self.work_manager.n_workers 

            #Slices this iteration into n_workers groups of segments, submits them to wm, splices results back together
            assignments, trajlabels = self.assign_iteration(nsegs[iiter], npts[iiter], parent_ids,
                                                            state_map, last_labels, pcoords, n_workers, n_iter)

            ##Do stuff with this iteration's results
                
            last_labels = trajlabels[:,-1].copy()
            assignments_ds[iiter, 0:nsegs[iiter], 0:npts[iiter]] = assignments
            if self.states:
                trajlabels_ds[iiter, 0:nsegs[iiter], 0:npts[iiter]]  = trajlabels
            
            del pcoords, assignments, parent_ids, trajlabels

        if sys.stdout.isatty() and not westpa.rc.quiet_mode:
            print('')


if __name__ == '__main__':
    WAssign().main()

