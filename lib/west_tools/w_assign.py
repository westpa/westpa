# Copyright (C) 2013 Matthew C. Zwier, Nick Rego, and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function, division; __metaclass__ = type
import sys
import logging
import math
from numpy import index_exp

from west.data_manager import seg_id_dtype, weight_dtype
from westpa.binning import index_dtype, assign_and_label, accumulate_labeled_populations 
from westtools import (WESTParallelTool, WESTDataReader, WESTDSSynthesizer, BinMappingComponent, 
                       ProgressIndicatorComponent)
import numpy
import westpa
from westpa import h5io
from westpa.h5io import WESTPAH5File
from westpa.extloader import get_object

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

def _assign_label_pop(n_iter, lb, ub, mapper, nstates, state_map, last_labels, parent_id_dsspec, weight_dsspec, pcoord_dsspec):

    nbins = len(state_map)-1
    parent_ids = parent_id_dsspec.get_iter_data(n_iter,index_exp[lb:ub])
    weights = weight_dsspec.get_iter_data(n_iter,index_exp[lb:ub])
    pcoords = pcoord_dsspec.get_iter_data(n_iter,index_exp[lb:ub])
    
    assignments, trajlabels, statelabels = assign_and_label(lb, ub, parent_ids,
                                               mapper.assign, nstates, state_map, last_labels, pcoords)
    pops = numpy.zeros((nstates+1,nbins+1), weight_dtype)
    accumulate_labeled_populations(weights, assignments, trajlabels, pops)
    return (assignments, trajlabels, pops, lb, ub, statelabels)

class WAssign(WESTParallelTool):
    prog='w_assign'
    description = '''\
Assign walkers to bins, producing a file (by default named "assign.h5")
which can be used in subsequent analysis.

For consistency in subsequent analysis operations, the entire dataset
must be assigned, even if only a subset of the data will be used. This
ensures that analyses that rely on tracing trajectories always know the
originating bin of each trajectory.


-----------------------------------------------------------------------------
Source data
-----------------------------------------------------------------------------

Source data is provided either by a user-specified function
(--construct-dataset) or a list of "data set specifications" (--dsspecs).
If neither is provided, the progress coordinate dataset ''pcoord'' is used.

To use a custom function to extract or calculate data whose probability
distribution will be calculated, specify the function in standard Python
MODULE.FUNCTION syntax as the argument to --construct-dataset. This function
will be called as function(n_iter,iter_group), where n_iter is the iteration
whose data are being considered and iter_group is the corresponding group
in the main WEST HDF5 file (west.h5). The function must return data which can
be indexed as [segment][timepoint][dimension].

To use a list of data set specifications, specify --dsspecs and then list the
desired datasets one-by-one (space-separated in most shells). These data set
specifications are formatted as NAME[,file=FILENAME,slice=SLICE], which will
use the dataset called NAME in the HDF5 file FILENAME (defaulting to the main
WEST HDF5 file west.h5), and slice it with the Python slice expression SLICE
(as in [0:2] to select the first two elements of the first axis of the
dataset). The ``slice`` option is most useful for selecting one column (or
more) from a multi-column dataset, such as arises when using a progress
coordinate of multiple dimensions.


-----------------------------------------------------------------------------
Specifying macrostates
-----------------------------------------------------------------------------

Optionally, kinetic macrostates may be defined in terms of sets of bins.
Each trajectory will be labeled with the kinetic macrostate it was most
recently in at each timepoint, for use in subsequent kinetic analysis.
This is required for all kinetics analysis (w_kintrace and w_kinmat).

There are three ways to specify macrostates:
  
  1. States corresponding to single bins may be identified on the command
     line using the --states option, which takes multiple arguments, one for
     each state (separated by spaces in most shells). Each state is specified
     as a coordinate tuple, with an optional label prepended, as in
     ``bound:1.0`` or ``unbound:(2.5,2.5)``. Unlabeled states are named
     ``stateN``, where N is the (zero-based) position in the list of states
     supplied to --states.
     
  2. States corresponding to multiple bins may use a YAML input file specified
     with --states-from-file. This file defines a list of states, each with a
     name and a list of coordinate tuples; bins containing these coordinates
     will be mapped to the containing state. For instance, the following
     file::

        ---
        states:
          - label: unbound
            coords:
              - [9.0, 1.0]
              - [9.0, 2.0]
          - label: bound
            coords:
              - [0.1, 0.0]

     produces two macrostates: the first state is called "unbound" and
     consists of bins containing the (2-dimensional) progress coordinate
     values (9.0, 1.0) and (9.0, 2.0); the second state is called "bound"
     and consists of the single bin containing the point (0.1, 0.0).
     
  3. Arbitrary state definitions may be supplied by a user-defined function,
     specified as --states-from-function=MODULE.FUNCTION. This function is
     called with the bin mapper as an argument (``function(mapper)``) and must
     return a list of dictionaries, one per state. Each dictionary must contain
     a vector of coordinate tuples with key "coords"; the bins into which each
     of these tuples falls define the state. An optional name for the state
     (with key "label") may also be provided.


-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "assign.h5") contains the following
attributes datasets:

  ``nbins`` attribute
    *(Integer)* Number of valid bins. Bin assignments range from 0 to
    *nbins*-1, inclusive.

  ``nstates`` attribute
    *(Integer)* Number of valid macrostates (may be zero if no such states are
    specified). Trajectory ensemble assignments range from 0 to *nstates*-1,
    inclusive, when states are defined.

  ``/assignments`` [iteration][segment][timepoint]
    *(Integer)* Per-segment and -timepoint assignments (bin indices).

  ``/npts`` [iteration]
    *(Integer)* Number of timepoints in each iteration.

  ``/nsegs`` [iteration]
    *(Integer)* Number of segments in each iteration.

  ``/labeled_populations`` [iterations][state][bin]
    *(Floating-point)* Per-iteration and -timepoint bin populations, labeled
    by most recently visited macrostate. The last state entry (*nstates-1*)
    corresponds to trajectories initiated outside of a defined macrostate.

  ``/bin_labels`` [bin]
    *(String)* Text labels of bins.

When macrostate assignments are given, the following additional datasets are
present:

  ``/trajlabels`` [iteration][segment][timepoint]
    *(Integer)* Per-segment and -timepoint trajectory labels, indicating the
    macrostate which each trajectory last visited.

  ``/state_labels`` [state]
    *(String)* Labels of states.

  ``/state_map`` [bin]
    *(Integer)* Mapping of bin index to the macrostate containing that bin.
    An entry will contain *nbins+1* if that bin does not fall into a 
    macrostate.
    
Datasets indexed by state and bin contain one more entry than the number of
valid states or bins. For *N* bins, axes indexed by bin are of size *N+1*, and
entry *N* (0-based indexing) corresponds to a walker outside of the defined bin
space (which will cause most mappers to raise an error). More importantly, for
*M* states (including the case *M=0* where no states are specified), axes
indexed by state are of size *M+1* and entry *M* refers to trajectories
initiated in a region not corresponding to a defined macrostate.

Thus, ``labeled_populations[:,:,:].sum(axis=1)[:,:-1]`` gives overall per-bin
populations, for all defined bins and 
``labeled_populations[:,:,:].sum(axis=2)[:,:-1]`` gives overall
per-trajectory-ensemble populations for all defined states. 

    
-----------------------------------------------------------------------------
Parallelization
-----------------------------------------------------------------------------

This tool supports parallelized binning, including reading/calculating input
data.


-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''
    
    def __init__(self):
        super(WAssign,self).__init__()
        
        # Parallel processing by default (this is not actually necessary, but it is
        # informative!)
        self.wm_env.default_work_manager = self.wm_env.default_parallel_work_manager
        
        self.data_reader = WESTDataReader()
        self.dssynth = WESTDSSynthesizer(default_dsname='pcoord')
        self.binning = BinMappingComponent()
        self.progress = ProgressIndicatorComponent()
        self.output_file = None
        self.output_filename = None
        self.states = []
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.binning.add_args(parser, suppress=['--bins-from-h5file'])
        self.dssynth.add_args(parser)
        
        sgroup = parser.add_argument_group('macrostate definitions').add_mutually_exclusive_group()
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

        agroup = parser.add_argument_group('other options')
        agroup.add_argument('-o', '--output', dest='output', default='assign.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)

        with self.data_reader:
            self.dssynth.h5filename = self.data_reader.we_h5filename
            self.dssynth.process_args(args)
            self.binning.process_args(args)

        if args.states:
            self.parse_cmdline_states(args.states)
        elif args.states_from_file:
            self.load_state_file(args.states_from_file)
        elif args.states_from_function:
            self.load_states_from_function(get_object(args.states_from_function,path=['.']))

        if self.states and len(self.states) < 2:
            raise ValueError('zero, two, or more macrostates are required')

        #self.output_file = WESTPAH5File(args.output, 'w', creating_program=True)
        self.output_filename = args.output
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


    def assign_iteration(self, n_iter, nstates, nbins, state_map, last_labels):
        ''' Method to encapsulate the segment slicing (into n_worker slices) and parallel job submission
            Submits job(s), waits on completion, splices them back together
            Returns: assignments, trajlabels, pops for this iteration'''

        futures = []

        iter_group = self.data_reader.get_iter_group(n_iter)
        nsegs, npts = iter_group['pcoord'].shape[:2]
        n_workers = self.work_manager.n_workers or 1
        assignments = numpy.empty((nsegs, npts), dtype=index_dtype)
        trajlabels = numpy.empty((nsegs, npts), dtype=index_dtype)
        statelabels = numpy.empty((nsegs, npts), dtype=index_dtype)
        pops = numpy.zeros((nstates+1,nbins+1), dtype=weight_dtype)

        #Submit jobs to work manager
        blocksize = nsegs // n_workers
        if nsegs % n_workers > 0:
            blocksize += 1

        def task_gen():
            if __debug__:
                checkset = set()
            for lb in xrange(0, nsegs, blocksize):
                ub = min(nsegs, lb+blocksize)
                if __debug__:
                    checkset.update(set(xrange(lb,ub)))
                args = ()
                kwargs = dict(n_iter=n_iter,
                              lb=lb, ub=ub, mapper=self.binning.mapper, nstates=nstates, state_map=state_map,
                              last_labels=last_labels, 
                              parent_id_dsspec=self.data_reader.parent_id_dsspec, 
                              weight_dsspec=self.data_reader.weight_dsspec,
                              pcoord_dsspec=self.dssynth.dsspec)
                yield (_assign_label_pop, args, kwargs)

                #futures.append(self.work_manager.submit(_assign_label_pop, 
                #kwargs=)
            if __debug__:
                assert checkset == set(xrange(nsegs)), 'segments missing: {}'.format(set(xrange(nsegs)) - checkset)

        #for future in self.work_manager.as_completed(futures):
        for future in self.work_manager.submit_as_completed(task_gen(), queue_size=self.max_queue_len):
            assign_slice, traj_slice, slice_pops, lb, ub, state_slice = future.get_result(discard=True)
            assignments[lb:ub, :] = assign_slice
            trajlabels[lb:ub, :] = traj_slice
            statelabels[lb:ub, :] = state_slice
            pops += slice_pops
            del assign_slice, traj_slice, slice_pops, state_slice

        del futures
        return (assignments, trajlabels, pops, statelabels)

    def go(self):
        assert self.data_reader.parent_id_dsspec._h5file is None
        assert self.data_reader.weight_dsspec._h5file is None
        if hasattr(self.dssynth.dsspec, '_h5file'):
            assert self.dssynth.dsspec._h5file is None
        pi = self.progress.indicator
        pi.operation = 'Initializing'
        with pi, self.data_reader, WESTPAH5File(self.output_filename, 'w', creating_program=True) as self.output_file:
            assign = self.binning.mapper.assign

            # We always assign the entire simulation, so that no trajectory appears to start
            # in a transition region that doesn't get initialized in one.
            iter_start = 1 
            iter_stop =  self.data_reader.current_iteration

            h5io.stamp_iter_range(self.output_file, iter_start, iter_stop)

            nbins = self.binning.mapper.nbins
            self.output_file.attrs['nbins'] = nbins 

            state_map = numpy.empty((self.binning.mapper.nbins+1,), index_dtype)
            state_map[:] = 0 # state_id == nstates => unknown state

            # Recursive mappers produce a generator rather than a list of labels
            # so consume the entire generator into a list
            labels = [label for label in self.binning.mapper.labels]

            self.output_file.create_dataset('bin_labels', data=labels, compression=9)

            if self.states:
                nstates = len(self.states)
                state_map[:] = nstates # state_id == nstates => unknown state
                state_labels = [state['label'] for state in self.states]

                for istate, sdict in enumerate(self.states):
                    assert state_labels[istate] == sdict['label'] #sanity check
                    state_assignments = assign(sdict['coords'])
                    for assignment in state_assignments:
                        state_map[assignment] = istate
                self.output_file.create_dataset('state_map', data=state_map, compression=9, shuffle=True)
                self.output_file['state_labels'] = state_labels #+ ['(unknown)']
            else:
                nstates = 0
            self.output_file.attrs['nstates'] = nstates

            iter_count = iter_stop - iter_start
            nsegs = numpy.empty((iter_count,), seg_id_dtype)
            npts = numpy.empty((iter_count,), seg_id_dtype)

            # scan for largest number of segments and largest number of points
            pi.new_operation ('Scanning for segment and point counts', iter_stop-iter_start)
            for iiter, n_iter in enumerate(xrange(iter_start,iter_stop)):
                iter_group = self.data_reader.get_iter_group(n_iter)
                nsegs[iiter], npts[iiter] = iter_group['pcoord'].shape[0:2]
                pi.progress += 1
                del iter_group

            pi.new_operation('Preparing output')

            # create datasets
            self.output_file.create_dataset('nsegs', data=nsegs, shuffle=True, compression=9)
            self.output_file.create_dataset('npts', data=npts, shuffle=True, compression=9)

            max_nsegs = nsegs.max()
            max_npts = npts.max()

            assignments_shape = (iter_count,max_nsegs,max_npts)
            assignments_dtype = numpy.min_scalar_type(nbins)
            assignments_ds = self.output_file.create_dataset('assignments', dtype=assignments_dtype, shape=assignments_shape,
                                                             compression=4, shuffle=True,
                                                             chunks=h5io.calc_chunksize(assignments_shape, assignments_dtype),
                                                             fillvalue=nbins)
            if self.states:
                trajlabel_dtype = numpy.min_scalar_type(nstates)
                trajlabels_ds = self.output_file.create_dataset('trajlabels', dtype=trajlabel_dtype, shape=assignments_shape,
                                                                compression=4, shuffle=True,
                                                                chunks=h5io.calc_chunksize(assignments_shape, trajlabel_dtype),
                                                                fillvalue=nstates)
                statelabels_ds = self.output_file.create_dataset('statelabels', dtype=trajlabel_dtype, shape=assignments_shape,
                                                                compression=4, shuffle=True,
                                                                chunks=h5io.calc_chunksize(assignments_shape, trajlabel_dtype),
                                                                fillvalue=nstates)

            pops_shape = (iter_count,nstates+1,nbins+1)
            pops_ds = self.output_file.create_dataset('labeled_populations', dtype=weight_dtype, shape=pops_shape,
                                                      compression=4, shuffle=True,
                                                      chunks=h5io.calc_chunksize(pops_shape, weight_dtype))
            h5io.label_axes(pops_ds, ['iteration', 'state', 'bin'])

            pi.new_operation('Assigning to bins', iter_stop-iter_start)
            last_labels = None # mapping of seg_id to last macrostate inhabited      
            for iiter, n_iter in enumerate(xrange(iter_start,iter_stop)):
                #get iteration info in this block

                if iiter == 0:
                    last_labels = numpy.empty((nsegs[iiter],), index_dtype)
                    last_labels[:] = nstates #unknown state

                #Slices this iteration into n_workers groups of segments, submits them to wm, splices results back together
                assignments, trajlabels, pops, statelabels = self.assign_iteration(n_iter, nstates, nbins, state_map, last_labels)

                ##Do stuff with this iteration's results

                last_labels = trajlabels[:,-1].copy()
                assignments_ds[iiter, 0:nsegs[iiter], 0:npts[iiter]] = assignments
                pops_ds[iiter] = pops
                if self.states:
                    trajlabels_ds[iiter, 0:nsegs[iiter], 0:npts[iiter]]  = trajlabels
                    statelabels_ds[iiter, 0:nsegs[iiter], 0:npts[iiter]]  = statelabels

                pi.progress += 1
                del assignments, trajlabels, pops, statelabels

            for dsname in 'assignments', 'npts', 'nsegs', 'labeled_populations', 'statelabels':
                h5io.stamp_iter_range(self.output_file[dsname], iter_start, iter_stop)

if __name__ == '__main__':
    WAssign().main()

