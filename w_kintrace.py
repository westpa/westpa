'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
import sys
import numpy
from westpa import h5io

import westpa
from west.data_manager import weight_dtype

from westpa.binning import index_dtype
from westpa.kinetics import find_macrostate_transitions
from westpa.kinetics._kinetics import _fast_transition_state_copy #@UnresolvedImport

ed_list_dtype = numpy.dtype([('istate', numpy.uint16), ('fstate', numpy.uint16), ('duration', numpy.float64),
                             ('weight', numpy.float64)])

log = logging.getLogger('westtools.w_kintrace')

class WKinTrace(WESTTool):
    prog='w_kintrace'
    description = '''\
Calculate state-to-state rates by tracing trajectories.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).
'''
    
    def __init__(self):
        super(WKinTrace,self).__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection() 
        self.output_file = None
        self.assignments_file = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
        
        parser.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')
        parser.add_argument('-o', '--output', dest='output', default='kintrace.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        
    def process_args(self, args):
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        self.data_reader.process_args(args)
        self.data_reader.open('r')
        self.iter_range.process_args(args)
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments do not span the requested iterations')
        
    def go(self): 
        nstates = self.assignments_file.attrs['nstates']
        start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter
                        
        durations_ds = self.output_file.create_dataset('durations', 
                                                       shape=(iter_count,0), maxshape=(iter_count,None),
                                                       dtype=ed_list_dtype,
                                                       chunks=(1,15360),
                                                       shuffle=True,
                                                       compression=9)
        durations_count_ds = self.output_file.create_dataset('duration_count',
                                                             shape=(iter_count,), dtype=numpy.int_, shuffle=True,compression=9)
        trace_fluxes_ds = self.output_file.create_dataset('trace_macro_fluxes', shape=(iter_count,nstates,nstates), dtype=weight_dtype,
                                                          chunks=h5io.calc_chunksize((iter_count,nstates,nstates),weight_dtype),
                                                          shuffle=True,
                                                          compression=9)

        
        # Put nice labels on things
        for ds in (self.output_file, durations_count_ds, trace_fluxes_ds):
            h5io.stamp_iter_range(ds, start_iter, stop_iter)
            
        # Calculate instantaneous rate matrices and trace trajectories
        last_state = None
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rIteration {}'.format(n_iter),end='')
                sys.stdout.flush()
            
            # Get data from the main HDF5 file
            iter_group = self.data_reader.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            nsegs, npts = iter_group['pcoord'].shape[0:2] 
            weights = seg_index['weight']
            parent_ids = seg_index['parent_id']
            
            # Get bin and traj. ensemble assignments from the previously-generated assignments file
            assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
            bin_assignments = numpy.require(self.assignments_file['assignments'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                            dtype=index_dtype)
            label_assignments = numpy.require(self.assignments_file['trajlabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                              dtype=index_dtype)
            
            # Prepare to run analysis
            macro_fluxes = numpy.zeros((nstates,nstates), weight_dtype)
            durations = []            

            # Estimate macrostate fluxes and calculate event durations using trajectory tracing
            # state is opaque to the find_macrostate_transitions function            
            state = _fast_transition_state_copy(iiter, nstates, parent_ids, last_state)
            find_macrostate_transitions(nstates, weights, label_assignments, 1.0/(npts-1), state, macro_fluxes, durations)
            last_state = state
            
            # Store trace-based kinetics data
            trace_fluxes_ds[iiter] = macro_fluxes
            durations_count_ds[iiter] = len(durations)
            if len(durations) > 0:
                durations_ds.resize((iter_count, max(len(durations), durations_ds.shape[1])))
                durations_ds[iiter,:len(durations)] = durations
                    
            # Do a little manual clean-up to prevent memory explosion
            del iter_group, weights, parent_ids, bin_assignments, label_assignments, state, macro_fluxes
        print()
        
if __name__ == '__main__':
    WKinTrace().main()
