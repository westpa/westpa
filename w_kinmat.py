'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
from collections import deque
import sys
import numpy
from westpa import h5io

import westpa
from west.data_manager import weight_dtype

from westpa.binning import index_dtype
from westpa.kinetics.matrates import estimate_rates

log = logging.getLogger('westtools.w_kinmat')

class WKinMat(WESTTool):
    prog='w_kinmat'
    description = '''\
Calculate populations, fluxes, and rates from weighted ensemble data
using a matrix approach to extrapolate to equilibrium values. This
analysis is only appropriate for simulations performed without sources
and sinks.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).
'''
    
    def __init__(self):
        super(WKinMat,self).__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection() 
        self.output_file = None
        self.assignments_file = None
        self.window_size = None
        self.all_lags = False
    
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
        parser.add_argument('--all-lags', action='store_true', default=False,
                            help='Use all possible lags within window of WINDOWSIZE')
        
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
        self.all_lags = bool(args.all_lags)
        
    #def calc_rate_matrices(self):
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

        
        for ds in (self.output_file, labeled_bin_pops_ds,labeled_bin_fluxes_ds):
            h5io.stamp_iter_range(ds, start_iter, stop_iter)
            
        h5io.label_axes(labeled_bin_pops_ds, ['iteration','state','bin'])
        for ds in (labeled_bin_fluxes_ds,):
            h5io.label_axes(ds, ['iteration','initial state','final state','inital bin','final bin'])

        # Calculate instantaneous rate matrices and trace trajectories
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
            weights_ring.append(weights)
            parent_ids_ring.append(parent_ids)
            bin_assignments_ring.append(bin_assignments)
            label_assignments_ring.append(label_assignments)
            
            # Estimate rates using bin-to-bin fluxes
            fluxes = estimate_rates(nbins, state_labels,
                                    weights_ring, parent_ids_ring, bin_assignments_ring, label_assignments_ring, state_map,
                                    self.all_lags)
            
            # Store bin-based kinetics data
            labeled_bin_fluxes_ds[iiter] = fluxes
                                
            # Do a little manual clean-up to prevent memory explosion
            del iter_group, weights, parent_ids, bin_assignments, label_assignments, fluxes
        print()
        
if __name__ == '__main__':
    WKinMat().main()
