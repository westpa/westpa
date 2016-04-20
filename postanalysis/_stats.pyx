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

# A cythoned version of the original function of the stats_process function,
# based on _kinetics.pyx

from __future__ import print_function,division
import cython
import numpy
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csgraph
import warnings
from collections import Counter
cimport numpy
cimport numpy as np
#cimport scipy.sparse as sp

ctypedef numpy.uint16_t index_t
ctypedef numpy.float64_t weight_t
ctypedef numpy.uint8_t bool_t
ctypedef numpy.int64_t trans_t
ctypedef numpy.uint_t uint_t # 32 bits on 32-bit systems, 64 bits on 64-bit systems

weight_dtype = numpy.float64  
index_dtype = numpy.uint16
bool_dtype = numpy.bool_

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef stats_process(numpy.ndarray[index_t, ndim=2] bin_assignments,
                    numpy.ndarray[weight_t, ndim=1] weights, 
                    numpy.ndarray[weight_t, ndim=2] fluxes, 
                    numpy.ndarray[weight_t, ndim=1] populations, 
                    numpy.ndarray[trans_t, ndim=2] trans, 
                    numpy.ndarray[index_t, ndim=2] mask,
                    str interval='timepoint'                        ):
    cdef:
        Py_ssize_t i,k
        index_t ibin,fbin,nsegs,npts
    nsegs = bin_assignments.shape[0]
    npts = bin_assignments.shape[1]

    if interval == 'timepoint':
        for i in xrange(0,npts - 1):
            for k in xrange(nsegs):
                ibin = bin_assignments[k,i]
                fbin = bin_assignments[k, i + 1]

                if mask[k, 0] == 1:
                    continue

                w = weights[k]

                fluxes[ibin, fbin] += w
                trans[ibin, fbin] += 1
                populations[ibin] += w
        return

    if interval == 'iteration':
        for k in xrange(nsegs):
            ibin = bin_assignments[k,i]
            fbin = bin_assignments[k, npts - 1]

            if mask[k, 0] == 1:
                continue

            w = weights[k]

            fluxes[ibin, fbin] += w
            trans[ibin, fbin] += 1
            populations[ibin] += w
        return

def normalize(m):
    nm = m.copy()

    row_sum = m.sum(1)
    ii = np.nonzero(row_sum)[0]
    nm[ii,:] = m[ii,:] / row_sum[ii][:, np.newaxis]

    return nm

def reweight_for_c(h5file, indices, nstates, nbins, state_labels, state_map, nfbins, istate, jstate, obs_threshold=1, total_fluxes=None, total_obs=None):

    # Instead of pulling in start and stop, we'll pull in a list of indices.
    # This way, it should support the bootstrap.
    total_fluxes, total_obs, total_pop = accumulate_statistics_list(h5file, indices, nfbins, total_fluxes, total_obs)

    flux_matrix = total_fluxes.copy()
    flux_matrix[total_obs < obs_threshold] = 0.0
    transition_matrix = normalize(flux_matrix)

    rw_bin_probs = steadystate_solve(transition_matrix)

    bin_last_state_map = np.tile(np.arange(nstates, dtype=np.int), nbins)
    bin_state_map = np.repeat(state_map[:-1], nstates)

    rw_color_probs = np.bincount(bin_last_state_map, weights=rw_bin_probs) 
    rw_state_probs = np.bincount(bin_state_map, weights=rw_bin_probs)

    rw_bin_transition_matrix = transition_matrix

    ii = np.nonzero(transition_matrix)

    rw_state_flux = calc_state_flux(rw_bin_transition_matrix[ii], ii[0], ii[1], rw_bin_probs, 
            bin_last_state_map, bin_state_map, nstates)

    # In truth, let's just return the rates for now.
    #return rw_state_flux, rw_color_probs, rw_state_probs, rw_bin_probs, rw_bin_transition_matrix
    return (rw_state_flux[istate,jstate] / (rw_color_probs[istate] / (rw_color_probs[istate] + rw_color_probs[jstate])))

def accumulate_statistics_list(h5file, iterations, nbins, total_fluxes=None, total_obs=None):
    if total_fluxes is None:
        assert total_obs is None
        total_fluxes = np.zeros((nbins, nbins), weight_dtype)
        total_obs = np.zeros((nbins, nbins), np.int64)

    rows = []
    cols = []
    obs = []
    flux = []
    total_pop = 0

    for iiter in iterations:
        iter_grp = h5file['iterations']['iter_{:08d}'.format(iiter)]

        rows.append(iter_grp['rows'][...])
        cols.append(iter_grp['cols'][...])
        obs.append(iter_grp['obs'][...])
        flux.append(iter_grp['flux'][...])
        total_pop += np.sum(h5file['bin_populations'][iiter, :], axis=0)

    rows, cols, obs, flux = map(np.hstack, [rows, cols, obs, flux])

    total_fluxes += sp.coo_matrix((flux, (rows, cols)), shape=(nbins, nbins)).todense()
    total_obs += sp.coo_matrix((obs, (rows, cols)), shape=(nbins, nbins)).todense()


    return total_fluxes, total_obs, total_pop

def calc_state_flux(trans_matrix, index1, index2, bin_probs, bin_last_state_map, bin_state_map, nstates):
    state_flux = np.zeros((nstates, nstates), np.float64)
    
    n_trans = index1.shape[0]
    for k in xrange(n_trans):
        ii = bin_last_state_map[index1[k]]
        jj = bin_state_map[index2[k]]

        if jj != nstates:
            state_flux[ii, jj] += trans_matrix[k] * bin_probs[index1[k]]

    return state_flux

def steadystate_solve(K):
    # Reformulate K to remove sink/source states
    n_components, component_assignments = csgraph.connected_components(K, connection="strong")
    largest_component = Counter(component_assignments).most_common(1)[0][0]
    components = np.where(component_assignments == largest_component)[0]

    ii = np.ix_(components, components)
    K_mod = K[ii]
    K_mod = normalize(K_mod)

    eigvals, eigvecs = np.linalg.eig(K_mod.T)
    eigvals = np.real(eigvals)
    eigvecs = np.real(eigvecs)

    maxi = np.argmax(eigvals)
    if not np.allclose(np.abs(eigvals[maxi]), 1.0):
        print('WARNING: Steady-state undetermined for current iteration')
        bin_prob = K.diagonal().copy()
        bin_prob = bin_prob / np.sum(bin_prob)
        return bin_prob

    sub_bin_prob = eigvecs[:, maxi] / np.sum(eigvecs[:, maxi])

    bin_prob = np.zeros(K.shape[0])
    bin_prob[components] = sub_bin_prob

    return bin_prob
