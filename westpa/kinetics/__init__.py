'''
Kinetics analysis library
'''

import logging
log = logging.getLogger(__name__)


from rate_averaging import RateAverager

import _kinetics
from _kinetics import (accumulate_labeled_populations, calculate_labeled_fluxes, #@UnresolvedImport
                       calculate_labeled_fluxes_alllags, #@UnresolvedImport
                       nested_to_flat_matrix, nested_to_flat_vector, #@UnresolvedImport
                       flat_to_nested_matrix, flat_to_nested_vector, find_macrostate_transitions) #@UnresolvedImport

"""
Internally, "labeled" objects (bin populations labeled by history, rate matrix elements labeled
by history) are stored as nested arrays -- e.g. rates[initial_label, final_label, initial_bin, final_bin].
These are converted to the flat forms required for, say, eigenvalue calculations internally, and the 
results converted back. This is because these conversions are not expensive, and saves users of 
this code from having to know how the flattened indexing works (something I screwed up all too
easily during development) -- mcz
"""

"""
Bin populations are treated separately from labeled bin populations because trajectories 
initiating in a transition region do not contribute to the labeled bin populations. -- mcz
"""

import numpy
import scipy.linalg
import scipy.sparse.linalg


from west.data_manager import weight_dtype
from collections import namedtuple

def get_steady_state_mult(T):
    T = numpy.asmatrix(T)
    for i in xrange(T.shape[0]):
        rowsum = T[i].sum()
        if rowsum > 0:
            T[i] /= rowsum
        else:
            T[i] = 0
    ivec = numpy.ones((T.shape[0],), weight_dtype) / T.shape[0]
    ivec.shape = (ivec.shape[0],1)
    ivec = numpy.asmatrix(ivec)
    ovec = ivec
    for i in xrange(50):
        T *= T
        ovec = T*ivec
    
    ss = numpy.array(ovec)
    ss.shape = (ss.shape[0],)
    ss /= ss.sum()
    return ss
        
        

def get_steady_state(rates):
    rates = rates.copy()#numpy.matrix(rates.copy())
    for i in xrange(rates.shape[0]):
        rowsum = rates[i].sum()
        if rowsum > 0:
            rates[i] = rates[i] / rowsum
        else:
            rates[i] = 0

    try:            
        vals, vecs = scipy.sparse.linalg.eigs(rates.T, k=1, which='LM')#,maxiter=100000,tol=1.0e-12)
        print(vals)
    except scipy.sparse.linalg.eigen.arpack.ArpackNoConvergence as e:
        return None
        
    ss = numpy.abs(vecs[:,0])
    #ss[ss<0]=0.0
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
        
    # Find steady-state solution
    ss = get_steady_state(rates)
    if ss is None:
        log.warning('no well-defined steady state; using average populations')
        ss = nested_to_flat_vector(labeled_pops)
    macro_rates = numpy.zeros((nstates,nstates), weight_dtype)
    
    # Sum over bins contributing to each state
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            for ibin in xrange(nbins):
                for jbin in xrange(nbins):
                    #sspop = ss[istate*nstates+ibin]
                    #rateelem = rates[istate*nbins+ibin,jstate*nbins+jbin]
                    sspop = ss[ibin*nstates+istate]
                    rateelem = rates[ibin*nstates+istate,jbin*nstates+jstate]
                    #print('state {}->{} bin {}->{} pop {} rate {}'.format(istate,jstate,ibin,jbin,sspop,rateelem))
                    macro_rates[istate,jstate] += sspop*rateelem
    
    # Normalize by total population in each trajectory ensemble
    for istate in xrange(nstates):
        traj_ens_pop = labeled_pops[istate].sum()
        macro_rates[istate] /= traj_ens_pop

    return flat_to_nested_vector(nstates, nbins, ss), macro_rates        


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

