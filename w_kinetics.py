'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
from collections import deque
from copy import deepcopy
import sys
import numpy, h5py
from numpy import index_exp
from westtools import h5io
import numba

import westpa
from west.data_manager import weight_dtype

from westpa.binning import UNKNOWN_INDEX
from westpa.kinetics import estimate_rates, nested_to_flat_matrix

ed_list_dtype = numpy.dtype([('istate', numpy.uint16), ('fstate', numpy.uint16), ('duration', numpy.float64),
                             ('weight', numpy.float64)])


log = logging.getLogger('westtools.w_kinetics')

@numba.jit(argtypes=[numba.int64, numba.int64, numba.float64, numba.float64, numba.float64, 
                     numba.uint16[:], numba.float64[:], numba.float64[:], numba.float64[:,:],
                     numba.float64[:,:], numba.object_],
           locals=dict(ipt=numba.int64,tm=numba.float64,flabel=numba.uint16,ilabel=numba.uint16,iistate=numba.uint16))
def _find_inner(nstates, npts, dt, weight, itime, label_assignments, last_entries, last_exits, last_completions, macro_fluxes, durations):
    # transitions never occur between the (overlapping) end point of previous iteration and beginning of
    # current iteration, so it suffices to start looking at timepoint 1 (and backwards to timepoint 0)
    for ipt in range(1,npts):
        tm = itime + ipt*dt
        flabel = label_assignments[ipt]
        ilabel = label_assignments[ipt-1]
        
        # if we have wound up in a new kinetic macrostate
        if flabel != ilabel:
            for iistate in xrange(nstates):
                if iistate != flabel:
                    macro_fluxes[iistate,flabel] += weight
                    t_ed = tm - last_exits[iistate]
                    durations.append((iistate,flabel,t_ed,weight))
                last_completions[iistate,flabel] = tm
            last_exits[ilabel] = tm
            last_entries[flabel] = tm
    return tm

def find_transitions(nstates, weights, label_assignments, dt, state, macro_fluxes, durations):
    # returns durations, new_state
    
    nsegs, npts = label_assignments.shape
    
    for seg_id in xrange(nsegs):
        state[seg_id]['last_time'] = _find_inner(nstates, npts, dt, weights[seg_id], 
                                                 state[seg_id]['last_time'], label_assignments[seg_id],
                                                 state[seg_id]['last_entries'], state[seg_id]['last_exits'],
                                                 state[seg_id]['last_completions'], macro_fluxes, durations)

class WKinetics(WESTTool):
    prog='w_kinetics'
    description = '''\
Calculate populations, fluxes, and rates from weighted ensemble data.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).
'''
    
    def __init__(self):
        super(WKinetics,self).__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection() 
        self.output_file = None
        self.assignments_file = None
        self.window_size = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
        
        parser.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')
        parser.add_argument('-o', '--output', dest='output', default='kinetics.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        parser.add_argument('-w', '--windowsize', type=int, default=1,
                            help='''Estimate kinetics over a maximum of WINDOWSIZE iterations.
                            (Default: %(default)s).''')

        
    def process_args(self, args):
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        self.data_reader.process_args(args)
        self.data_reader.open('r')
        self.iter_range.process_args(args)
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.window_size = args.windowsize
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments do not span the requested iterations')
        
    #def calc_rate_matrices(self):
    @profile
    def go(self):
        nbins = self.assignments_file.attrs['nbins']
        state_labels = self.assignments_file['state_labels'][...]
        state_map = self.assignments_file['state_map'][...]
        nstates = len(state_labels)
        start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter
        
        weights_ring = deque(maxlen=self.window_size)
        parent_ids_ring = deque(maxlen=self.window_size)
        bin_assignments_ring = deque(maxlen=self.window_size)
        label_assignments_ring = deque(maxlen=self.window_size)
        
        labeled_vector_shape = (iter_count,nstates,nbins)
        labeled_matrix_shape = (iter_count,nstates,nstates,nbins,nbins)
        
        bin_pops_ds = self.output_file.create_dataset('bin_pops', shape=(iter_count,nbins), dtype=weight_dtype)
        labeled_bin_pops_ds = self.output_file.create_dataset('labeled_bin_pops',
                                                              shape=labeled_vector_shape,
                                                              chunks=h5io.calc_chunksize(labeled_vector_shape, weight_dtype),
                                                              compression=9,
                                                              dtype=weight_dtype)
        labeled_bin_fluxes_ds = self.output_file.create_dataset('labeled_bin_fluxes',
                                                                shape=labeled_matrix_shape,
                                                                chunks=h5io.calc_chunksize(labeled_matrix_shape, weight_dtype),
                                                                compression=9,
                                                                dtype=weight_dtype)
        labeled_bin_rates_ds = self.output_file.create_dataset('labeled_bin_rates',
                                                               shape=labeled_matrix_shape,
                                                               chunks=h5io.calc_chunksize(labeled_matrix_shape, weight_dtype),
                                                               compression=9,
                                                               dtype=weight_dtype)

        avg_macro_fluxes = numpy.zeros((nstates,nstates), weight_dtype)
        
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

        
        for ds in (self.output_file, bin_pops_ds,labeled_bin_pops_ds,labeled_bin_fluxes_ds,labeled_bin_rates_ds,
                   durations_count_ds, trace_fluxes_ds):
            h5io.stamp_iter_range(ds, start_iter, stop_iter)
            
        h5io.label_axes(bin_pops_ds, ['iteration', 'bin'])
        h5io.label_axes(labeled_bin_pops_ds, ['iteration','state','bin'])
        for ds in (labeled_bin_fluxes_ds,labeled_bin_rates_ds):
            h5io.label_axes(ds, ['iteration','initial state','final state','inital bin','final bin'])

        # Calculate instantaneous rate matrices and trace trajectories
        last_state = None
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rIteration {}'.format(n_iter),end='')
                sys.stdout.flush()
            
            iter_group = self.data_reader.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            nsegs, npts = iter_group['pcoord'].shape[0:2] 
            weights = seg_index['weight']
            parent_ids = seg_index['parent_id']
            
            assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
            
            bin_assignments = self.assignments_file['assignments'][assignment_iiter + numpy.s_[:nsegs,:npts]]
            label_assignments = self.assignments_file['trajlabels'][assignment_iiter + numpy.s_[:nsegs,:npts]]
            
            macro_fluxes = numpy.zeros((nstates,nstates), weight_dtype)
            durations = []            
            weights_ring.append(weights)
            parent_ids_ring.append(parent_ids)
            bin_assignments_ring.append(bin_assignments)
            label_assignments_ring.append(label_assignments)
            
            # Estimate rates using bin-to-bin fluxes
            rateest = estimate_rates(nbins, state_labels,
                                     weights_ring, parent_ids_ring, bin_assignments_ring, label_assignments_ring, state_map)

            # Estimate macrostate fluxes and calculate event durations using trajectory tracing

            # this needs some "fast state copy" function -- 70% of run time between here and "end" comment
            state = numpy.empty((nsegs,), numpy.object_)            
            for seg_id, parent_id in enumerate(parent_ids):
                if iiter == 0 or parent_id < 0:
                    state[seg_id] = dict(last_entries = numpy.zeros((nstates,), numpy.float64),
                                         last_exits = numpy.zeros((nstates,), numpy.float64),
                                         last_completions = numpy.zeros((nstates,nstates), numpy.float64),
                                         last_time = 0.0)
                else:
                    parent_state = last_state[parent_id]
                    state[seg_id] = dict(last_time = parent_state['last_time'], last_entries = parent_state['last_entries'].copy(), last_exits = parent_state['last_exits'].copy(),last_completions = parent_state['last_completions'].copy())
            del last_state
            # end "fast state copy"
                        
            find_transitions(nstates, weights, label_assignments, 1.0, state, macro_fluxes, durations)
            last_state = state
            
            durations_count_ds[iiter] = len(durations)
            if len(durations) > 0:
                durations_ds.resize((iter_count, max(len(durations), durations_ds.shape[1])))
                durations_ds[iiter,:len(durations)] = durations
                
            total_pops = rateest.labeled_bin_pops.sum(axis=1)
            for istate in xrange(nstates):
                macro_fluxes[istate,:] /= total_pops[istate]
            
            trace_fluxes_ds[iiter] = macro_fluxes
            avg_macro_fluxes += macro_fluxes

            bin_pops_ds[iiter] = rateest.bin_pops
            labeled_bin_pops_ds[iiter] = rateest.labeled_bin_pops
            labeled_bin_fluxes_ds[iiter] = rateest.labeled_bin_fluxes
            labeled_bin_rates_ds[iiter] = rateest.labeled_bin_rates
            print(macro_fluxes)
        
            del iter_group, weights, parent_ids, bin_assignments, label_assignments, rateest, state, macro_fluxes
        print()
        
if __name__ == '__main__':
    WKinetics().main()
