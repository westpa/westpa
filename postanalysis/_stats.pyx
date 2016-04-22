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
ctypedef unsigned short Ushort

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
#cpdef weight_t[:,:] normalize(weight_t[:,:] m):
cpdef normalize(numpy.ndarray[weight_t, ndim=2] tm):

    cdef:
        weight_t row_sum
        Py_ssize_t x, y, nfbins
        numpy.ndarray[weight_t, ndim=2] m

    # We should send this in, rather than calculating it here.
    # We also want to convert to memoryviews, ultimately, if possible.
    # However, doing that means we pass nothing BUT memoryviews, which I'm not yet ready for.
    m = tm.copy()
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
        int[:] _rows, _cols, _obs, _ins, _nrows, _ncols, _nobs
        weight_t[:] _flux, _total_pop, rw_bin_probs, _nflux
        int n_trans, _nstates
        long[:,:] _ii

        weight_t[:,:] _total_fluxes
        int[:,:] _total_obs



    # CREATE NUMPY ARRAYS
    # This is a temporary measure that fixes some segfaults, which implies I'm probably off by
    # a little bit.  Memory heavy, but whatever.
    # It breaks depending on things, so I need to root that out.  Clearly, nnz is larger than that.
    nnz = len(rows)*2
    total_fluxes = np.zeros((nfbins, nfbins), weight_dtype)
    total_obs = np.zeros((nfbins, nfbins), intc_dtype)
    transition_matrix = np.zeros((nfbins, nfbins), weight_dtype)
    nrows = np.zeros((nnz), intc_dtype)
    ncols = np.zeros((nnz), intc_dtype)
    nobs = np.zeros((nnz), intc_dtype)
    nflux = np.zeros((nnz), weight_dtype)


    # CREATE MEMORYVIEWS
    # These are for what we sent in...
    _rows = rows
    _cols = cols
    _obs = obs
    _flux = flux
    _ins = insert
    _nstates = nstates
    # ... these are for functions we'll be using.
    _nrows = nrows
    _ncols = ncols
    _nobs = nobs
    _nflux = nflux
    _total_fluxes = total_fluxes
    _total_obs = total_obs


    # We want to reconstruct the missing datasets, so we assume we have the first part of the block and
    # reconstruct appropriately.
    # There is probably a better way to do this.
    new_indices = np.zeros(((indices.shape[0]*stride)), dtype=indices.dtype)
    indices = regenerate_subsampled_indices(indices, new_indices, len(indices), stride)


    _total_fluxes, _total_obs = accumulate_statistics_list(_rows, _cols, _obs, _flux, _ins, indices, nfbins, nnz, total_fluxes, total_obs, nrows, ncols, nobs, nflux, len(indices))

    # This is still an issue.  Clean up the handling of numpy arrays vs. memoryviews.
    total_fluxes[total_obs < obs_threshold] = 0.0
    transition_matrix = np.asarray(total_fluxes.copy())
    transition_matrix = np.asarray(normalize(transition_matrix))

    rw_bin_probs = np.asarray(steadystate_solve(transition_matrix.copy()))

    bin_last_state_map = np.tile(np.arange(nstates, dtype=np.int), nbins)
    bin_state_map = np.repeat(state_map[:-1], nstates)

    rw_color_probs = np.bincount(bin_last_state_map, weights=rw_bin_probs) 
    rw_state_probs = np.bincount(bin_state_map, weights=rw_bin_probs)

    rw_bin_transition_matrix = transition_matrix

    ii = np.nonzero(transition_matrix)

    n_trans = ii[0].shape[0]
    state_flux = np.zeros((nstates, nstates), np.float64)
    rw_state_flux = calc_state_flux(rw_bin_transition_matrix[ii], ii[0], ii[1], rw_bin_probs, 
            bin_last_state_map, bin_state_map, nstates, state_flux, n_trans)

    if rw_color_probs[istate] != 0.0:
        return (rw_state_flux[istate,jstate] / (rw_color_probs[istate] / (rw_color_probs[istate] + rw_color_probs[jstate])))
    else:
        return 0.0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef regenerate_subsampled_indices(Ushort[:] iin, Ushort[:] iout, int ilen, int stride):

    cdef:
        int i, si

    with nogil:
        for i in range(ilen):
            for si in range(stride):
                iout[(i*stride)+si] = iin[i] + si

    return iout

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef accumulate_statistics_list(int[:] hrows, int[:] hcols, int[:] hobs, weight_t[:] hflux, int[:] hins, index_t[:] iterations, Py_ssize_t nbins, Py_ssize_t nnz, weight_t[:,:] total_fluxes, int[:,:] total_obs, int[:] rows, int[:] cols, int[:] obs, weight_t[:] flux, int itermax):

    cdef:
        index_t curriter, elem, iiter, ilem, ipop

    curriter = 0

    with nogil:
        for iter in range(itermax):
            iiter = iterations[iter]
            for ilem in range(hins[iiter], hins[iiter+1]):
                rows[curriter] = hrows[ilem]
                cols[curriter] = hcols[ilem]
                obs[curriter] = hobs[ilem]
                flux[curriter] = hflux[ilem]
                curriter += 1

        total_fluxes = dense_to_sparse_float(rows, cols, flux, curriter, total_fluxes)
        total_obs = dense_to_sparse_int(rows, cols, obs, curriter, total_obs)


    return total_fluxes, total_obs

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef weight_t[:,:] dense_to_sparse_float(int[:] rows, int[:] cols, weight_t[:] flux, int elem, weight_t[:,:] total_fluxes) nogil:

    cdef:
        int i

    for i in range(elem):
        total_fluxes[rows[i], cols[i]] += flux[i]

    return total_fluxes

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int[:,:] dense_to_sparse_int(int[:] rows, int[:] cols, int[:] obs, int elem, int[:,:] total_obs) nogil:

    cdef:
        int i

    for i in range(elem):
        total_obs[rows[i], cols[i]] += obs[i]

    return total_obs



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef weight_t[:,:] calc_state_flux(weight_t[:] trans_matrix, long[:] index1, long[:] index2, weight_t[:] bin_probs, long[:] bin_last_state_map, Ushort[:] bin_state_map, int nstates, weight_t[:,:] state_flux, int n_trans) nogil:
    
    cdef:
        int k, ii, jj


    for k in xrange(n_trans):
        ii = bin_last_state_map[index1[k]]
        jj = bin_state_map[index2[k]]

        if jj != nstates:
            state_flux[ii, jj] += (trans_matrix[k] * bin_probs[index1[k]])

    return state_flux

def steadystate_solve(K):

    cdef:
        weight_t[:] _bin_prob
        weight_t[:,:] K_mod

    bin_prob = np.zeros(K.shape[0])
    _bin_prob = bin_prob
    # Reformulate K to remove sink/source states
    n_components, component_assignments = csgraph.connected_components(K, connection="strong")
    largest_component = Counter(component_assignments).most_common(1)[0][0]
    components = np.where(component_assignments == largest_component)[0]

    ii = np.ix_(components, components)
    K_mod = K[ii].copy()
    K_mod = normalize(np.asarray(K_mod))

    eigvals, eigvecs = np.linalg.eig(K_mod.T)
    eigvals = np.real(eigvals)
    eigvecs = np.real(eigvecs)

    maxi = np.argmax(eigvals)
    if not np.allclose(np.abs(eigvals[maxi]), 1.0):
        bin_prob = K.diagonal().copy()
        bin_prob = bin_prob / np.sum(bin_prob)
        return bin_prob

    sub_bin_prob = eigvecs[:, maxi] / np.sum(eigvecs[:, maxi])

    bin_prob[components] = sub_bin_prob

    return bin_prob
