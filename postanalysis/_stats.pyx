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
from scipy.sparse import csgraph
import warnings
from collections import Counter
cimport numpy
cimport numpy as np
cimport scipy.linalg.cython_lapack as cl
cimport scipy.linalg
import scipy.linalg

ctypedef numpy.uint16_t index_t
ctypedef numpy.float64_t weight_t
ctypedef numpy.uint8_t bool_t
ctypedef numpy.int64_t trans_t
ctypedef numpy.uint_t uint_t # 32 bits on 32-bit systems, 64 bits on 64-bit systems
ctypedef unsigned short Ushort
ctypedef double complex Cdouble

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
cpdef int normalize(weight_t[:,:] m, Py_ssize_t nfbins) nogil:

    cdef:
        weight_t row_sum
        Py_ssize_t x, y

    for y in range(nfbins):
        row_sum = 0
        for x in range(nfbins):
            row_sum += m[y,x]
        if row_sum != 0:
            for x in range(nfbins):
                m[y,x] /= row_sum
    return 0


#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.cdivision(True)
def reweight_for_c(rows, cols, obs, flux, insert, indices, nstates, nbins, state_labels, state_map, nfbins, istate, jstate, stride, bin_last_state_map, bin_state_map, obs_threshold=1):



    # Instead of pulling in start and stop, we'll pull in a list of indices.
    # This way, it should support the bootstrap.
    cdef:
        int[:] _rows, _cols, _obs, _ins, _nrows, _ncols, _nobs
        weight_t[:] _flux, _total_pop, _rw_bin_probs, _nflux
        int n_trans, _nstates, lind, _nfbins, _stride, _obs_threshold, nnz, nlind
        long[:,:] _ii
        Ushort[:] _indices
        Ushort[:] _new_indices

        weight_t[:,:] _total_fluxes, _transition_matrix, _rw_state_flux
        int[:,:] _total_obs
        #double[:] eigvals, eigvalsi
        #double[:,:] eigvecs, WORK



    # CREATE NUMPY ARRAYS
    # This is a temporary measure that fixes some segfaults, which implies I'm probably off by
    # a little bit.  Memory heavy, but whatever.
    # It breaks depending on things, so I need to root that out.  Clearly, nnz is larger than that.
    _flux = flux
    nnz = len(flux)
    lind = indices.shape[0]
    nlind = indices.shape[0]*stride
    total_fluxes = np.zeros((nfbins, nfbins), weight_dtype)
    total_obs = np.zeros((nfbins, nfbins), intc_dtype)
    transition_matrix = np.zeros((nfbins, nfbins), weight_dtype)
    strong_transition_matrix = np.zeros((nfbins, nfbins), weight_dtype)
    rw_bin_probs = np.zeros(nfbins, weight_dtype)
    new_indices = np.zeros(((nlind)), dtype=indices.dtype)
    rw_state_flux = np.zeros((nstates, nstates), np.float64)
    state_flux = np.zeros((nstates, nstates), np.float64)
    eigvals = np.zeros((nfbins), np.float64)
    eigvalsi = np.zeros((nfbins), np.float64)
    eigvecs = np.zeros((nfbins, nfbins), np.float64)
    WORK = np.zeros((nfbins*4, nfbins*4), np.float64)


    # CREATE MEMORYVIEWS
    # These are for what we sent in...
    _rows = rows
    _cols = cols
    _obs = obs
    _flux = flux
    _ins = insert
    _nstates = nstates
    # ... these are for functions we'll be using.
    _total_fluxes = total_fluxes
    _total_obs = total_obs
    _transition_matrix = transition_matrix
    _nfbins = nfbins
    _stride = stride
    _obs_threshold = obs_threshold
    _indices = indices
    _new_indices = new_indices
    _rw_state_flux = rw_state_flux
    _rw_bin_probs = rw_bin_probs


    #NOGIL
    # Reconstruct dataset.  We're just passing the same thing back and forth between functions.
    #with nogil:
    regenerate_subsampled_indices(indices, new_indices, lind, stride)
    #accumulate_statistics_list(_rows, _cols, _obs, _flux, _ins, new_indices, _nfbins, nnz, transition_matrix, total_obs, nlind)
    #_transition_matrix, _total_obs = accumulate_statistics_list(rows, cols, obs, flux, insert, new_indices, nfbins, nnz, transition_matrix, total_obs, nlind)
    accumulate_fluxes(rows, cols, obs, flux, insert, new_indices, nnz, transition_matrix, nlind)
    accumulate_obs(rows, cols, obs, flux, insert, new_indices, nnz, total_obs, nlind)

    remove_under_obs(transition_matrix, total_obs, obs_threshold, nfbins)
    normalize(_transition_matrix, _nfbins)

    #print("new round!")
    steadystate_solve(transition_matrix, strong_transition_matrix, rw_bin_probs, nfbins, eigvals, eigvalsi, eigvecs, WORK)
    #print(np.real(eigvals))
    #print(np.real(eigvecs))

    # Can kill this calculation by sending in the data already...
    #bin_last_state_map = np.tile(np.arange(nstates, dtype=np.int), nbins)
    #bin_state_map = np.repeat(state_map[:-1], nstates)

    rw_color_probs = np.bincount(bin_last_state_map, weights=rw_bin_probs) 
    rw_state_probs = np.bincount(bin_state_map, weights=rw_bin_probs)


    ii = np.nonzero(transition_matrix)

    n_trans = ii[0].shape[0]
    calc_state_flux(transition_matrix[ii], ii[0], ii[1], rw_bin_probs, 
            bin_last_state_map, bin_state_map, nstates, rw_state_flux, n_trans)

    #print(rw_color_probs)
    #print(rw_state_flux)
    if rw_color_probs[istate] != 0.0:
        #print(rw_state_flux[istate,jstate] / (rw_color_probs[istate] / (rw_color_probs[istate] + rw_color_probs[jstate])))
        return (rw_state_flux[istate,jstate] / (rw_color_probs[istate] / (rw_color_probs[istate] + rw_color_probs[jstate])))
    else:
        return 0.0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int regenerate_subsampled_indices(Ushort[:] iin, Ushort[:] iout, int ilen, int stride) nogil:

    cdef:
        int i, si

    for i in range(ilen):
        for si in range(stride):
            iout[(i*stride)+si] = iin[i] + si

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int accumulate_fluxes(int[:] hrows, int[:] hcols, int[:] hobs, weight_t[:] hflux, int[:] hins, index_t[:] iterations, Py_ssize_t nnz, weight_t[:,:] total_fluxes, int itermax) nogil:

    cdef:
        index_t curriter, elem, iiter, ilem, ipop

    curriter = 0

    for iter in range(itermax):
        iiter = iterations[iter]
        for ilem in range(hins[iiter], hins[iiter+1]):
            total_fluxes[hrows[ilem], hcols[ilem]] += hflux[ilem]



    # We're just operating on them in memory, right?
    #return total_fluxes, total_obs
    #return total_fluxes
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int accumulate_obs(int[:] hrows, int[:] hcols, int[:] hobs, weight_t[:] hflux, int[:] hins, index_t[:] iterations, Py_ssize_t nnz, int[:,:] total_obs, int itermax) nogil:

    cdef:
        index_t curriter, elem, iiter, ilem, ipop

    curriter = 0

    for iter in range(itermax):
        iiter = iterations[iter]
        for ilem in range(hins[iiter], hins[iiter+1]):
            if ilem < nnz and iiter+1 < itermax:
                total_obs[hrows[ilem], hcols[ilem]] += hobs[ilem]



    # We're just operating on them in memory, right?
    #return total_fluxes, total_obs
    #return total_obs
    return 0

#@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int dense_to_sparse_float(int[:] rows, int[:] cols, weight_t[:] flux, int elem, weight_t[:,:] total_fluxes):

    cdef:
        int i

    for i in range(elem):
        total_fluxes[rows[i], cols[i]] += flux[i]

    return 0

#@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int dense_to_sparse_int(int[:] rows, int[:] cols, int[:] obs, int elem, int[:,:] total_obs):

    cdef:
        int i

    for i in range(elem):
        total_obs[rows[i], cols[i]] += obs[i]

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int remove_under_obs(weight_t[:,:] flux, int[:,:] obs, int threshold, int nbins) nogil:

    cdef:
        int x, y

    for x in range(nbins):
        for y in range(nbins):
            if obs[x,y] < threshold:
                flux[x,y] = 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int calc_state_flux(weight_t[:] trans_matrix, long[:] index1, long[:] index2, weight_t[:] bin_probs, long[:] bin_last_state_map, Ushort[:] bin_state_map, int nstates, weight_t[:,:] state_flux, int n_trans) nogil:
    
    cdef:
        int k, ii, jj


    for k in xrange(n_trans):
        ii = bin_last_state_map[index1[k]]
        jj = bin_state_map[index2[k]]

        if jj != nstates:
            state_flux[ii, jj] += (trans_matrix[k] * bin_probs[index1[k]])

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int steadystate_solve(weight_t[:,:] K, weight_t[:,:] K_mod, weight_t[:] bin_prob, int K_shape, double[:] eigvals, double[:] eigvalsi, double[:,:] eigvecs, double[:,:] WORK):

    cdef:
        #double[:] eigvals, eigvalsi
        #double[:,:] eigvecs, WORK
        #double[:] eigvals_o
        #double[:,:] eigvecs_o
        long[:] components
        double max, eigsum
        int n_components, components_assignments, largest_component, maxi, x, y, n, INFO, LWORK

        # POINTERS

        int  *_INFO, *_K_shape, *_LWORK
        double *_K_mod, *_eigvals, *_eigvecs, *_WORK, *_eigvalsi, *_eigvals_o, *_eigvalsi_o
    _K_shape = &K_shape
    _INFO = &INFO
    _LWORK = &LWORK
    INFO = 0
    LWORK = K_shape * 4
    _K_mod = &K_mod[0,0]
    _eigvals = &eigvals[0]
    #_eigvals_o = &eigvals_o[0]
    #_eigvalsi_0 = &eigvalsi_o[0]
    _eigvecs = &eigvecs[0,0]
    _eigvalsi = &eigvalsi[0]
    _WORK = &WORK[0,0]

    # CREATE NUMPY STRUCTURES
    # Reformulate K to remove sink/source states
    n_components, component_assignments = csgraph.connected_components(K, connection="strong")
    largest_component = Counter(component_assignments).most_common(1)[0][0]
    components = np.where(component_assignments == largest_component)[0]
    n = len(components)
    #maxi = 0

    # This works.
    for x in range(n):
        for y in range(n):
            K_mod[components[x], components[y]] = K[components[x], components[y]]
    normalize(K_mod, K_shape)

    #eigvals, eigvecs = np.linalg.eig(K_mod.T)
    #eigvals, eigvecs = cl.dgeev(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    # JOBVL: 'V' left eigenvectors
    # JOBVR: 'N' right eigenvectors
    # N: Matrix order.  K_shape
    # A: Transition matrix
    # LDA: Leading dimension of the array.  Also K_shape, as we're square.
    # WR, WI: Real and imaginary eigenvalues.
    # VL: Left eigenvectors
    # LDVL: leading dimension of the array VL.  Should be K_shape
    # VR: Right eigenvectors
    # LDVR: same as for LDVL
    # WORK: 
    eigsum = 0.0

    with nogil:
        cl.dgeev('N', 'V', _K_shape, _K_mod, _K_shape, _eigvals, _eigvalsi, _eigvecs, _K_shape, _eigvecs, _K_shape, _WORK, _LWORK, _INFO)
        # Second _eigvecs can probably be null.
        for x in range(K_shape):
            if x == 0:
                max = eigvals[0]
                maxi = x
            else:
                if max < eigvals[x]:
                    max = eigvals[x]
                    maxi = x
        #if not np.allclose(np.abs(eigvals[maxi]), 1.0):
            # Need to fix.
            #bin_prob = K.diagonal().copy()
            #bin_prob = bin_prob / np.sum(bin_prob)
        #    return 0
            #return bin_prob

        # Apparently, for a memoryview object, you have to explicitly set everything by hand.  Which is adorable.
        # This works, as well.
        #for i in range(n):
        #    eigsum += eigvecs[components[i], maxi]
        #for i in range(n):
        #    bin_prob[components[i]] = eigvecs[components[i], maxi]
        #    bin_prob[components[i]] /= eigsum
        for i in range(n):
            eigsum += eigvecs[maxi, components[i]]
        for i in range(n):
            bin_prob[components[i]] = eigvecs[maxi, components[i]]
            bin_prob[components[i]] /= eigsum

    return 0
