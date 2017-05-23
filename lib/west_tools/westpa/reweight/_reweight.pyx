# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
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
from libc.math cimport isnan

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
            # Should this be 0?
            # .... this should super be 0.  What?
            ibin = bin_assignments[k,0]
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef reweight_for_c(rows, cols, obs, flux, insert, indices, nstates, nbins, state_labels, state_map, nfbins, istate, jstate, stride, bin_last_state_map, bin_state_map, return_obs, obs_threshold=1):



    # Instead of pulling in start and stop, we'll pull in a list of indices.
    # This way, it should support the bootstrap.
    cdef:
        int[:] _rows, _cols, _obs, _ins, _nrows, _ncols, _nobs, _visited
        long[:]  _bin_last_state_map
        weight_t[:] _flux, _total_pop, _rw_bin_probs, _nflux, _rw_color_probs, _rw_state_probs
        double[:] _eigvals, _eigvalsi
        int n_trans, _nstates, lind, _nfbins, _stride, _obs_threshold, nnz, nlind, i, j, _istate, _jstate
        Ushort[:] _indices, _bin_state_map
        weight_t _return_value

        Ushort[:] _new_indices

        weight_t[:,:] _total_fluxes, _transition_matrix, _rw_state_flux, _strong_transition_matrix
        double[:,:] _WORK, _eigvecs
        int[:,:] _total_obs, _graph
        #bint _return_flux, _return_states, _return_color
        str _return_obs
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
    visited = np.zeros((nfbins), intc_dtype)
    graph = np.zeros((nfbins, nfbins+1), dtype=intc_dtype)
    rw_bin_probs = np.zeros(nfbins, weight_dtype)
    new_indices = np.zeros(((nlind)), dtype=indices.dtype)
    rw_state_flux = np.zeros((nstates, nstates), np.float64)
    state_flux = np.zeros((nstates, nstates), weight_dtype)
    eigvals = np.zeros((nfbins), np.float64)
    eigvalsi = np.zeros((nfbins), np.float64)
    eigvecs = np.zeros((nfbins, nfbins), np.float64)
    WORK = np.zeros((nfbins*4, nfbins*4), np.float64)
    rw_color_probs = np.zeros((nstates), weight_dtype)
    rw_state_probs = np.zeros((nbins), weight_dtype)


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
    _strong_transition_matrix = strong_transition_matrix
    _nfbins = nfbins
    _stride = stride
    _obs_threshold = obs_threshold
    _indices = indices
    _new_indices = new_indices
    _rw_state_flux = rw_state_flux
    _rw_bin_probs = rw_bin_probs
    _eigvals = eigvals
    _eigvalsi = eigvalsi
    _eigvals = eigvals
    _eigvecs = eigvecs
    _WORK = WORK
    _graph = graph
    _visited = visited
    _rw_color_probs = rw_color_probs
    _rw_state_probs = rw_state_probs
    _bin_last_state_map = bin_last_state_map
    _bin_state_map = bin_state_map
    _istate = istate
    _jstate = jstate
    _return_obs = return_obs


    #NOGIL
    # Reconstruct dataset.  We're just passing the same thing back and forth between functions.
    with nogil:
        for i in range(_nfbins):
            for j in range(1, _nfbins+1):
                _graph[i, j] = _nfbins
        regenerate_subsampled_indices(_indices, _new_indices, lind, _stride)
        accumulate_fluxes(_rows, _cols, _obs, _flux, _ins, _new_indices, nnz, _transition_matrix, nlind)
        accumulate_obs(_rows, _cols, _obs, _flux, _ins, _new_indices, nnz, _total_obs, nlind)

        remove_under_obs(_transition_matrix, _total_obs, _obs_threshold, _nfbins)
        normalize(_transition_matrix, _nfbins)
        steadystate_solve(_transition_matrix, _strong_transition_matrix, _rw_bin_probs, _nfbins, _eigvals, _eigvalsi, _eigvecs, _WORK, _graph, _visited)

        for i in range(_nfbins):
            _rw_color_probs[_bin_last_state_map[i]] += _rw_bin_probs[i]
            _rw_state_probs[_bin_state_map[i]] += _rw_bin_probs[i]


        calc_state_flux(_transition_matrix, _rw_bin_probs, _bin_last_state_map, _bin_state_map, _nstates, _rw_state_flux, _nfbins)

    # This allows us to use the same function for all three types.
    # Return conditional fluxes.
    if _return_obs == b'F':
        _return_value = _rw_state_flux[_istate,_jstate]
        if isnan(_return_value) == True:
            return 0.0
        else:
            return _return_value
    # Return state probabilities.
    elif _return_obs == b'S':
        _return_value = _rw_state_probs[_istate]
        if isnan(_return_value) == True:
            return 0.0
        else:
            return _return_value
    # Return color (ensemble) probabilities
    elif _return_obs == b'C':
        _return_value = _rw_color_probs[_istate]
        if isnan(_return_value) == True:
            return 0.0
        else:
            return _return_value
    # Return the rates.
    elif _return_obs == b'R':
        if _rw_color_probs[_istate] != 0.0:
            _return_value = (_rw_state_flux[_istate,_jstate] / (_rw_color_probs[_istate] / (_rw_color_probs[_istate] + _rw_color_probs[_jstate])))
            if isnan(_return_value) == True:
                return 0.0
            else:
                return _return_value
        else:
            # We have no ensemble probability, and as such, cannot have a flux.
            return 0.0
    # Return the populations.
    elif _return_obs == b'P':
        return rw_bin_probs

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int regenerate_subsampled_indices(Ushort[:] iin, Ushort[:] iout, int ilen, int stride) nogil:

    cdef:
        int i, si

    # go over the range of all indices within iin
    for i in range(ilen):
        # Run over the length of the stride.
        for si in range(stride):
            iout[(i*stride)+si] = iin[i] + si

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int accumulate_fluxes(int[:] hrows, int[:] hcols, int[:] hobs, weight_t[:] hflux, int[:] hins, Ushort[:] iterations, Py_ssize_t nnz, weight_t[:,:] total_fluxes, int itermax) nogil:

    cdef:
        index_t curriter, elem, iiter, ipop
        long ilem

    curriter = 0

    for iter in range(itermax):
        iiter = iterations[iter]
        for ilem in range(hins[iiter], hins[iiter+1]):
            # Not sure if this is necessary, here...
            if ilem < nnz and iiter+1 < itermax:
                total_fluxes[hrows[ilem], hcols[ilem]] += hflux[ilem]

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int accumulate_obs(int[:] hrows, int[:] hcols, int[:] hobs, weight_t[:] hflux, int[:] hins, Ushort[:] iterations, Py_ssize_t nnz, int[:,:] total_obs, int itermax) nogil:

    cdef:
        index_t curriter, elem, iiter, ipop
        long ilem

    curriter = 0

    for iter in range(itermax):
        iiter = iterations[iter]
        for ilem in range(hins[iiter], hins[iiter+1]):
            if ilem < nnz and iiter+1 < itermax:
                total_obs[hrows[ilem], hcols[ilem]] += hobs[ilem]

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
cpdef int calc_state_flux(weight_t[:, :] trans_matrix, weight_t[:] bin_probs, long[:] bin_last_state_map, Ushort[:] bin_state_map, int nstates, weight_t[:,:] state_flux, int K_shape) nogil:
    
    cdef:
        int i, j, ii, jj


    for i in range(K_shape):
        for j in range(K_shape):

            ii = bin_last_state_map[i]
            jj = bin_state_map[j]

            if jj != nstates:
                state_flux[ii, jj] += (trans_matrix[i, j] * bin_probs[i])

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int steadystate_solve(weight_t[:,:] K, weight_t[:,:] K_mod, weight_t[:] bin_prob, int K_shape, double[:] eigvals, double[:] eigvalsi, double[:,:] eigvecs, double[:,:] WORK, int[:,:] graph, int[:] visited) nogil:

    cdef:
        int[:] components
        int[:,:] _graph
        double max, eigsum
        int n_components, components_assignments, largest_component, maxi, x, y, n, INFO, LWORK, i, j

        # POINTERS

        int  *_INFO, *_K_shape, *_LWORK
        double *_K_mod, *_eigvals, *_eigvecs, *_WORK, *_eigvalsi
    _K_shape = &K_shape
    _INFO = &INFO
    _LWORK = &LWORK
    INFO = 0
    LWORK = K_shape * 4
    _K_mod = &K_mod[0,0]
    _eigvals = &eigvals[0]
    _eigvecs = &eigvecs[0,0]
    _eigvalsi = &eigvalsi[0]
    _WORK = &WORK[0,0]
    _graph = graph

    for i in range(K_shape):
        if visited[i] == 0:
            visited[i] = 1
            return_strong_component(K, K_shape, _graph, i, i, visited)
    n = 0
    for i in range(K_shape):
        if graph[i, 0] >= graph[n, 0]:
            n = i
    # I suspect this may be giving us issues?
    #components = _graph[n, :K_shape+1]

    maxi = 0
    eigsum = 0.0
    # This all works!
    for x in range(K_shape):
        #i = components[x+1]
        i = graph[n, x+1]
        for y in range(K_shape):
            #j = components[y+1]
            j = graph[n, y+1]
            if i != K_shape and j != K_shape:
                K_mod[i, j] = K[i, j]
    normalize(K_mod, K_shape)
    cl.dgeev('N', 'V', _K_shape, _K_mod, _K_shape, _eigvals, _eigvalsi, _eigvecs, _K_shape, _eigvecs, _K_shape, _WORK, _LWORK, _INFO)
    for x in range(K_shape):
        if x == 0:
            max = eigvals[0]
            maxi = x
        else:
            if max < eigvals[x]:
                max = eigvals[x]
                maxi = x
    # We need to go over the whole range and pick out non K_shape elements.
    # This probably no longer needs to be done, now...
    for i in range(K_shape):
        #x = components[i+1]
        x = graph[n, i+1]
        if x != K_shape:
            #eigsum += eigvecs[maxi, components[i+1]]
            eigsum += eigvecs[maxi, x]
    for i in range(K_shape):
        #x = components[i+1]
        x = graph[n, i+1]
        if x != K_shape:
            #bin_prob[components[i+1]] = eigvecs[maxi, components[i+1]]
            #bin_prob[components[i+1]] /= eigsum
            bin_prob[x] = eigvecs[maxi, x]
            bin_prob[x] /= eigsum

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int return_strong_component(weight_t[:,:] K, int K_shape, int[:, :] graph, int i, int z, int[:] visited) nogil:

    cdef:
        int j, y 

    if graph[z, 0] == 0:
        graph[z, 0] += 1
        graph[z, 1] = i
    for j in xrange(K_shape):
        if i != j:
            if K[i, j] > 0.0:
                # Strongly connected!
                if visited[j] == 0:
                    graph[z, 0] += 1
                    y = graph[z, 0]
                    graph[z, y] = j
                    # We only want to call it when we haven't visited it before.
                    # We don't want to call, THEN modify and check.  Otherwise, we could be doing many calls.
                    visited[j] = 1
                    return_strong_component(K, K_shape, graph, j, z, visited)


    return 0
