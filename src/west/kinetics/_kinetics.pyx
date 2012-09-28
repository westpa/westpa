from __future__ import print_function,division
import cython
import numpy
cimport numpy

ctypedef numpy.uint16_t index_t
ctypedef numpy.float64_t weight_t
ctypedef numpy.uint8_t bool_t

weight_dtype = numpy.float64  
index_dtype = numpy.uint16
bool_dtype = numpy.bool_


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef flux_assign(numpy.ndarray[weight_t, ndim=1] weights,
                  numpy.ndarray[index_t, ndim=1] init_assignments,
                  numpy.ndarray[index_t, ndim=1] final_assignments,
                  numpy.ndarray[weight_t, ndim=2] flux_matrix):
    cdef:
        Py_ssize_t m,n
        index_t i, j
    n = len(weights)
    for m from 0 <= m < n:
        i = init_assignments[m]
        j = final_assignments[m]
        flux_matrix[i,j] += weights[m]
        
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef pop_assign(numpy.ndarray[weight_t, ndim=1] weights,
                 numpy.ndarray[index_t, ndim=1] assignments,
                 numpy.ndarray[weight_t, ndim=1] populations):
    cdef:
        Py_ssize_t m,n
        index_t i,
    n = len(weights)
    for m from 0 <= m < n:
        i = assignments[m]
        populations[i] += weights[m]

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef calc_rates(numpy.ndarray[weight_t, ndim=3] fluxes,
                 numpy.ndarray[weight_t, ndim=2] populations,
                 numpy.ndarray[weight_t, ndim=3] rates,
                 numpy.ndarray[bool_t,ndim=2,cast=True] masks):
    '''Calculate a series of rate matrices from a series of flux and population matrices.
    A set of boolean matrices, of the same shape as populations, is also produced, to be 
    used for generating masks for the rate matrices where initial state populations are
    zero.'''
    
    cdef:
        Py_ssize_t narrays, nbins
        index_t iarray, i, j
        
    narrays = fluxes.shape[0]
    nbins = fluxes.shape[1]
    
    for iarray from 0 <= iarray < narrays:
        for i from 0 <= i < nbins:
            if populations[iarray,i] == 0.0:
                masks[iarray,i] = 0
                for j from 0 <= j < nbins:
                    rates[iarray,i,j] = 0
            else:
                masks[iarray,i] = 1
                for j from 0 <= j < nbins:
                    rates[iarray,i,j] = fluxes[iarray,i,j] / populations[iarray,i]
    

