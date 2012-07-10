from __future__ import division
import numpy, sys, math, os
cimport numpy, cython

#from posix.unistd cimport size_t#
#from cpython.mem cimport PyMem_Malloc, PyMem_Free
#include "unistd.pxd"
#include "cpython/mem.pxd"
from libc.stdlib cimport RAND_MAX, rand, srand
include "posix/unistd.pxd"
include "cpython/mem.pxd"

ctypedef fused _fptype:
    numpy.float32_t
    numpy.float64_t 

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef _fptype _autocorrel_elem(numpy.ndarray[_fptype,ndim=1] xs, long k):
    '''Calculate the ``k``-lag sample autocorrelation of ``xs``'''
    cdef long N = xs.shape[0]
    cdef long i
    cdef double cdif
    cdef double xbar = xs.mean()
    cdef double norm = 0.0
    cdef double rho = 0.0

    with nogil:    
        for i in range(N):
            cdif = xs[i] - xbar
            norm += cdif*cdif
            
            if i < N-k:
                rho += cdif * (xs[i+k]-xbar)
     
        rho *= N / ((N-k)*norm)
    
    return <_fptype> rho
        
def autocorrel_elem(xs, k):
    if xs.dtype == numpy.float64:
        return _autocorrel_elem[numpy.float64_t](xs,k)
    elif xs.dtype == numpy.float32:
        return _autocorrel_elem[numpy.float32_t](xs,k)
    else:
        raise TypeError('unsupported type')


