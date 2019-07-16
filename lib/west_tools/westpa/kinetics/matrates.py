# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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

'''
Routines for implementing Letteri et al.'s macrostate-to-macrostate rate calculations
using extrapolation to steady-state populations from average rate matrices
'''

"""
Internally, "labeled" objects (bin populations labeled by history, rate matrix elements labeled
by history) are stored as nested arrays -- e.g. rates[initial_label, final_label, initial_bin, final_bin].
These are converted to the flat forms required for, say, eigenvalue calculations internally, and the 
results converted back. This is because these conversions are not expensive, and saves users of 
this code from having to know how the flattened indexing works (something I screwed up all too
easily during development) -- mcz
"""


import logging, warnings
log = logging.getLogger(__name__)

from _kinetics import (calculate_labeled_fluxes, #@UnresolvedImport
                       calculate_labeled_fluxes_alllags, #@UnresolvedImport
                       labeled_flux_to_rate, #@UnresolvedImport
                       nested_to_flat_matrix, nested_to_flat_vector, #@UnresolvedImport
                       flat_to_nested_vector, _reduce_labeled_rate_matrix_to_macro) #@UnresolvedImport

import numpy
import scipy.linalg

from west.data_manager import weight_dtype

class ConsistencyWarning(UserWarning):
    pass

def get_steady_state(rates):
    '''Get steady state solution for a rate matrix. As an optimization, returns the
    flattened labeled population vector (of length nstates*nbins); to convert to the
    nested vector used for storage, use nested_to_flat_vector().'''

    rates = rates.copy()

    # Convert to a transition probability matrix
    for i in range(rates.shape[0]):
        rowsum = rates[i,:].sum()
        if rowsum > 0:
            rates[i,:] = rates[i,:] / rowsum
        else:
            if rates[:,i].sum() != 0:
                warnings.warn('sink microstate in rate matrix', ConsistencyWarning)
            rates[i,:] = 0

    try:
        vals, vecs = scipy.linalg.eig(rates.T)
    except Exception:
        log.debug('exception obtaining eigenvectors', exc_info=True)
        return None

    vals = numpy.abs(vals)
    log.debug('eigenvalues: {!r}'.format(list(reversed(sorted(vals)))))
    asort = numpy.argsort(vals)
    vec = vecs[:,asort[-1]]
    ss = numpy.abs(vec)

    ss /= ss.sum()
    return ss

def get_macrostate_rates(labeled_rates, labeled_pops, extrapolate=True):
    '''Using a labeled rate matrix and labeled bin populations, calculate the steady state
    probability distribution and consequent state-to-state rates.

    Returns ``(ss, macro_rates)``, where ``ss`` is the steady-state probability distribution
    and ``macro_rates`` is the state-to-state rate matrix.'''

    nstates, nbins = labeled_pops.shape

    rates = nested_to_flat_matrix(labeled_rates)

    # Find steady-state solution
    if extrapolate:
        ss = get_steady_state(rates)
        if ss is None:
            warnings.warn('no well-defined steady state; using average populations',ConsistencyWarning)
            ss = nested_to_flat_vector(labeled_pops)
    else:
        ss = nested_to_flat_vector(labeled_pops)

    macro_rates = _reduce_labeled_rate_matrix_to_macro(nstates, nbins, rates, ss)

    return flat_to_nested_vector(nstates, nbins, ss), macro_rates

def estimate_rates(nbins, state_labels, weights, parent_ids, bin_assignments, label_assignments, state_map, labeled_pops,
                   all_lags=False,
                   labeled_fluxes = None, labeled_rates = None, unlabeled_rates = None):
    '''Estimate fluxes and rates over multiple iterations. The number of iterations is determined by how many
    vectors of weights, parent IDs, bin assignments, and label assignments are passed.

    If ``all_lags`` is true, then the average is over all possible lags within the length-N window given, otherwise
    simply the length N lag.

    Returns labeled flux matrix, labeled rate matrix, and unlabeled rate matrix.'''

    assert len(weights) == len(parent_ids) == len(bin_assignments) == len(label_assignments)
    nstates = len(state_labels)
    nbins = labeled_pops.shape[1]-1

    # Prepare output arrays
    if labeled_fluxes is None:
        labeled_fluxes = numpy.zeros((nstates, nstates, nbins, nbins), weight_dtype)
    else:
        labeled_fluxes.fill(0.0)

    if labeled_rates is None:
        labeled_rates = numpy.zeros_like(labeled_fluxes)
    else:
        labeled_rates.fill(0.0)
    
    if unlabeled_rates is None:
        unlabeled_rates = numpy.zeros((nbins,nbins), weight_dtype)
    else:
        unlabeled_rates.fill(0.0)


    # Loop over all possible windows to accumulate flux matrix
    # flux matrix is [initial_label][final_label][initial_bin][final_bin]
    if all_lags:
        twindow = calculate_labeled_fluxes_alllags(nstates, weights, parent_ids, bin_assignments, label_assignments, labeled_fluxes)
    else:
        twindow = calculate_labeled_fluxes(nstates, weights, parent_ids, bin_assignments, label_assignments, labeled_fluxes)
    labeled_fluxes /= twindow

    # Calculate rate matrix for this window, using populations from the last iteration (which correspond
    # to the weights that contribute to the flux matrix)
    labeled_flux_to_rate(labeled_fluxes, labeled_pops, labeled_rates)

    # Calculate an unlabeled rate matrix
    unlabeled_fluxes = numpy.sum(labeled_fluxes,axis=(0,1))
    unlabeled_pops = labeled_pops.sum(axis=0)
    unlabeled_rates[...] = labeled_flux_to_rate(unlabeled_fluxes[None,None,:,:], unlabeled_pops[None,:])[0,0]

    return labeled_fluxes, labeled_rates, unlabeled_rates

