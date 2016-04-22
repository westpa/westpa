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
import h5py
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
intc_dtype = numpy.intc

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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef normalize(numpy.ndarray[weight_t, ndim=2] tm):

    cdef:
        weight_t row_sum
        Py_ssize_t x, y, nfbins
        numpy.ndarray[weight_t, ndim=2] m

    m = tm.copy()
    # We should send this in, rather than calculating it here.
    # We also want to convert to memoryviews, ultimately, if possible.
    nfbins = m.shape[0]

    with nogil:
        for y in range(nfbins):
            row_sum = 0
            for x in range(nfbins):
                row_sum += m[y,x]
            if row_sum != 0:
                for x in range(nfbins):
                    m[y,x] /= row_sum

    return m

#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.cdivision(True)
def reweight_for_c(rows, cols, obs, flux, insert, indices, nstates, nbins, state_labels, state_map, nfbins, istate, jstate, stride, obs_threshold=1):



    # Instead of pulling in start and stop, we'll pull in a list of indices.
    # This way, it should support the bootstrap.
    cdef:
        int[:] _rows, _cols, _obs, _int
        weight_t[:] _flux, _total_pop, rw_bin_probs
    _rows = rows
    _cols = cols
    _obs = obs
    _flux = flux
    _ins = insert

    # We want to reconstruct the missing datasets, so we assume we have the first part of the block and
    # reconstruct appropriately.
    if stride > 1:
    #    print("Show me!")
        new_indices = np.repeat(indices, stride)
        new_indices = np.ones_like(new_indices)
        curriter = 0
        for ind in range(len(indices)):
            i = indices[ind]
            for x in range(stride):
                if curriter < len(indices)*stride:
                    if i+x < len(insert)-1:
                        new_indices[curriter] = i+x
                        curriter +=1
        indices = new_indices


    # This is a temporary measure that fixes some segfaults, which implies I'm probably off by
    # a little bit.  Memory heavy, but whatever.
    nnz = len(rows)*2
    indices = np.sort(indices)
    total_fluxes, total_obs = accumulate_statistics_list(_rows, _cols, _obs, _flux, _ins, indices, nfbins, nnz)

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
    # Truth be told, this shouldn't really happen.  But we should make sure we're not gonna blow out with NaNs.
    #if rw_color_probs[istate] == 0:
    #    return 0
    #else:
        #print(rw_state_flux[istate,jstate] / rw_color_probs[istate])
    #    return (rw_state_flux[istate,jstate] / (rw_color_probs[istate]))
    #return rw_state_flux[istate,jstate]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef accumulate_statistics_list(int[:] hrows, int[:] hcols, int[:] hobs, weight_t[:] hflux, int[:] hins, index_t[:] iterations, Py_ssize_t nbins, Py_ssize_t nnz):

    cdef:
        int[:] _rows, _cols, _obs, _ins
        weight_t[:] _flux
        index_t curriter, elem, iiter, ilem, ipop
        weight_t[:,:] _total_fluxes
        int[:,:] _total_obs

    # We need to know how many nonzero elements we have.
    # Bit nasty, since we'll need to iterate through and get them all in one place.
    # Ugh, basically.  Probably best to feed this in from the main function call.
    itermax = len(iterations)

    total_fluxes = np.zeros((nbins, nbins), weight_dtype)
    total_obs = np.zeros((nbins, nbins), intc_dtype)
    rows = np.zeros((nnz), intc_dtype)
    cols = np.zeros((nnz), intc_dtype)
    obs = np.zeros((nnz), intc_dtype)
    flux = np.zeros((nnz), weight_dtype)
    _rows = rows
    _cols = cols
    _obs = obs
    _flux = flux


    curriter = 0

    # Memory view for numpy arrays.
    _total_fluxes = total_fluxes
    _total_obs = total_obs

    with nogil:
        for iter in range(itermax):
            iiter = iterations[iter]
            for ilem in range(hins[iiter], hins[iiter+1]):
                #print(nnz, iter, iiter, ilem, curriter, iterations)
                _rows[curriter] = hrows[ilem]
                _cols[curriter] = hcols[ilem]
                _obs[curriter] = hobs[ilem]
                _flux[curriter] = hflux[ilem]
                curriter += 1
            # Sum along the axis, here...
            # Wait, we're not even using this,

    total_fluxes = sp.coo_matrix((_flux, (_rows, _cols)), shape=(nbins, nbins)).todense()
    total_obs = sp.coo_matrix((_obs, (_rows, _cols)), shape=(nbins, nbins)).todense()


    return total_fluxes, total_obs

def calc_state_flux(trans_matrix, index1, index2, bin_probs, bin_last_state_map, bin_state_map, nstates):
    
    cdef:
        weight_t[:,:] _trans_matrix, _state_flux
        weight_t[:] _bin_probs
        Py_ssize_t k

    state_flux = np.zeros((nstates, nstates), np.float64)
    _state_flux = state_flux
    _trans_matrix = trans_matrix
    _bin_probs = bin_probs

    n_trans = index1.shape[0]
    for k in xrange(n_trans):
        ii = bin_last_state_map[index1[k]]
        jj = bin_state_map[index2[k]]

        if jj != nstates:
            # This isn't working because we're not slicing the trans_matrix correctly.
            # Unsure (doubtful) that it's working correctly, currently.

            #_state_flux[ii, jj] += (_trans_matrix[0, k] * _bin_probs[index1[k]])
            #print(trans_matrix[0,k], bin_probs[index1[k]])
            _state_flux[ii, jj] += (trans_matrix[0, k] * bin_probs[index1[k]])

    return state_flux

def steadystate_solve(K):

    cdef:
        weight_t[:] _bin_prob

    bin_prob = np.zeros(K.shape[0])
    _bin_prob = bin_prob
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
        diag = K.diagonal().copy()
        for i in range(diag.shape[1]):
            _bin_prob[i] = diag[0, i]
        bin_prob = bin_prob / np.sum(bin_prob)
        # THIS PART CAN GET BENT ASSHOLE
        return bin_prob

    sub_bin_prob = eigvecs[:, maxi] / np.sum(eigvecs[:, maxi])

    bin_prob[components] = sub_bin_prob

    return bin_prob
