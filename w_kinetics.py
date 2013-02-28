'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
from collections import namedtuple, deque
import sys
import numpy, h5py
import scipy.linalg
import scipy.sparse.linalg
from westtools import h5io

import westpa
from westpa.binning.assign import index_dtype, UNKNOWN_INDEX
from westpa.kinetics._kinetics import (accumulate_labeled_populations, calculate_labeled_fluxes, #@UnresolvedImport
                                       calculate_labeled_fluxes_alllags, #@UnresolvedImport
                                       nested_to_flat_matrix, nested_to_flat_vector, #@UnresolvedImport
                                       flat_to_nested_matrix, flat_to_nested_vector) #@UnresolvedImport
from west.data_manager import seg_id_dtype, n_iter_dtype, weight_dtype

from work_managers import make_work_manager

import numba
from numba import f8, i8, u8 #@UnresolvedImport

log = logging.getLogger('westtools.w_kinetics')

"""
Internally, "labeled" objects (bin populations labeled by history, rate matrix elements labeled
by history) are stored as nested arrays -- e.g. rates[initial_label, final_label, initial_bin, final_bin].
These are converted to the flat forms required for, say, eigenvalue calculations internally, and the 
results converted back. This is because these conversions are not expensive, and saves users of 
this code from having to know how the flattened indexing works (something I screwed up all too
easily during development) -- mcz
"""

"""
Bin populations are stored separately from labeled bin populations because trajectories 
initiating in a transition region do not contribute to the labeled bin populations. -- mcz
"""

def get_steady_state(rates):
    rates = numpy.matrix(rates.copy())
    for i in xrange(rates.shape[0]):
        rowsum = rates[i].sum()
        if rowsum > 0:
            rates[i] = rates[i] / rowsum

    try:
        _, vecs = scipy.sparse.linalg.eigs(rates.T, k=1, which='LM')
    except scipy.sparse.linalg.eigen.arpack.ArpackError:
        # ill-defined steady-state
        ss = numpy.empty((rates.shape[0],), weight_dtype)
        ss[:] = float('nan')
    else:
        ss = numpy.abs(vecs[:,0])
        ss /= ss.sum()
    
    return ss

def labeled_flux_to_rate(labeled_fluxes, labeled_pops):
    '''Convert a labeled flux matrix and corresponding labeled bin populations to
    a labeled rate matrix.'''
            
    rates = numpy.empty_like(labeled_fluxes)
    for istate in xrange(labeled_fluxes.shape[0]):
        for ibin in xrange(labeled_fluxes.shape[1]):
            pop = labeled_pops[istate,ibin]
            if pop > 0:
                rates[istate,ibin] = labeled_fluxes[istate,ibin] / pop
            else:
                rates[istate,ibin] = 0
    
    return rates

def get_macrostate_rates(labeled_rates, labeled_pops):
    '''Using a labeled rate matrix and labeled bin populations, calculate the steady state
    probability distribution and consequent state-to-state rates.
    
    Returns ``(ss, macro_rates)``, where ``ss`` is the steady-state probability distribution
    and ``macro_rates`` is the state-to-state rate matrix.'''
    
    nstates, nbins = labeled_pops.shape
        
    rates = nested_to_flat_matrix(labeled_rates)
    flatpops = nested_to_flat_vector(labeled_pops)
    
    for i in xrange(rates.shape[0]):
        if flatpops[i] == 0:
            rates[i] = 0
        else:
            rates[i] /= flatpops[i]
    
    # Find steady-state solution
    ss = get_steady_state(rates)
    macro_rates = numpy.zeros((nstates,nstates), weight_dtype)
    
    # Sum over bins contributing to each state
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            for ibin in xrange(nbins):
                for jbin in xrange(nbins):
                    sspop = ss[istate*nstates+ibin]
                    rateelem = rates[istate*nbins+ibin,jstate*nbins+jbin]
                                                              
                    macro_rates[istate,jstate] += sspop*rateelem
    
    # Normalize by total population in each trajectory ensemble
    for istate in xrange(nstates):
        traj_ens_pop = labeled_pops[istate].sum()
        macro_rates[istate] /= traj_ens_pop

    return ss, macro_rates        


_rate_result = namedtuple('_rate_result', ('bin_pops','labeled_bin_pops','labeled_bin_fluxes', 'labeled_bin_rates'))

def estimate_rates(nbins, state_labels, weights, parent_ids, bin_assignments, label_assignments, state_map,
                   all_lags=False):
    '''Estimate rates over multiple iterations.
    Returns unlabeled and labeled bin populations, labeled flux matrix, and
    (the instantaneous estimate of) the rate matrix.'''
    
    assert len(weights) == len(parent_ids) == len(bin_assignments) == len(label_assignments)
    
    niters = len(weights)
    nstates = len(state_labels)
        
    # Estimate microstate and trajectory populations over entire window
    bin_pops = numpy.zeros((nbins,), weight_dtype)
    labeled_bin_pops = numpy.zeros((nstates,nbins), weight_dtype)
    
    for iiter in xrange(niters):
        accumulate_labeled_populations(weights[iiter], bin_assignments[iiter], label_assignments[iiter], bin_pops,
                                       labeled_bin_pops)
    
    bin_pops /= niters
    labeled_bin_pops /= niters 
        
    # Loop over all possible windows to accumulate flux matrix
    # flux matrix is [initial_label][final_label][initial_bin][final_bin]
    fluxes = numpy.zeros((nstates, nstates, nbins, nbins), weight_dtype)
    if all_lags:
        twindow = calculate_labeled_fluxes_alllags(nstates, weights, parent_ids, bin_assignments, label_assignments, fluxes)
    else:
        twindow = calculate_labeled_fluxes(nstates, weights, parent_ids, bin_assignments, label_assignments, fluxes)
    fluxes /= twindow
    
    rates = labeled_flux_to_rate(fluxes, labeled_bin_pops)    
    return _rate_result(bin_pops, labeled_bin_pops, fluxes, rates)

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
        self.output_file = None
        self.assignments_file = None
        self.window_size = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        
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
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.window_size = args.windowsize
        
    def go(self):
        nbins = self.assignments_file.attrs['nbins']
        state_labels = self.assignments_file['state_labels'][...]
        state_map = self.assignments_file['state_map'][...]
        nstates = len(state_labels)
        start_iter, stop_iter = h5io.get_iter_range(self.assignments_file)
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
        
        for ds in (bin_pops_ds,labeled_bin_pops_ds,labeled_bin_fluxes_ds,labeled_bin_rates_ds):
            h5io.stamp_iter_range(ds, start_iter, stop_iter)
            
        h5io.label_axes(bin_pops_ds, ['iteration', 'bin'])
        h5io.label_axes(labeled_bin_pops_ds, ['iteration','state','bin'])
        for ds in (labeled_bin_fluxes_ds,labeled_bin_rates_ds):
            h5io.label_axes(ds, ['iteration','initial state','final state','inital bin','final bin'])

        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rIteration {}'.format(n_iter),end='')
                sys.stdout.flush()
            
            iter_group = self.data_reader.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            nsegs, npts = iter_group['pcoord'].shape[0:2] 
            weights = seg_index['weight']
            parent_ids = seg_index['parent_id']
            bin_assignments = self.assignments_file['assignments'][iiter,:nsegs,:npts]
            label_assignments = self.assignments_file['trajlabels'][iiter,:nsegs,:npts]
            
            weights_ring.append(weights)
            parent_ids_ring.append(parent_ids)
            bin_assignments_ring.append(bin_assignments)
            label_assignments_ring.append(label_assignments)
            
            rateest = estimate_rates(nbins, state_labels,
                                     weights_ring, parent_ids_ring, bin_assignments_ring, label_assignments_ring, state_map)
            
            bin_pops_ds[iiter] = rateest.bin_pops
            labeled_bin_pops_ds[iiter] = rateest.labeled_bin_pops
            labeled_bin_fluxes_ds[iiter] = rateest.labeled_bin_fluxes
            labeled_bin_rates_ds[iiter] = rateest.labeled_bin_rates
        
            del iter_group, weights, parent_ids, bin_assignments, label_assignments, rateest
        print()
        

if __name__ == '__main__':
    WKinetics().main()
