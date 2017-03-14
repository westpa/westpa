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

from __future__ import division
import numpy, sys, math, os
cimport numpy, cython
from cpython.buffer cimport *
from numpy cimport *

from libc.math cimport floor, ceil, log10
#from libc.stdlib cimport RAND_MAX, rand, srand

ctypedef fused _fptype:
    numpy.float32_t
    numpy.float64_t 
    
# Number of bytes in a dataset above which correlation time estimates are 
# not performed by generating synthetic data sets once and then evaluating 
# the autocorrelation elements on them, but rather by generating a new
# synthetic dataset each time an autocorrelation element is required.
# Essentially, this limits the RAM used during the correlation time
# calculation to approximately CORRELTIME_CROSSOVER.  Keep in mind that
# the space savings is at the cost of needing *many* more random numbers,
# which will almost certainly become the rate-limiting step of the
# calculation.
cdef Py_ssize_t CORRELTIME_CROSSOVER = 512*1024*1024
    

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef _fptype _autocorrel_elem(_fptype[:] xs, long k):
    cdef:
        long N = xs.shape[0]
        long i
        double cdif
        double xbar = 0.0
        double norm = 0.0
        double rho = 0.0

    with nogil:
        # calculate the mean manually, because it's faster than calling out to
        # numpy.mean
        for i in range(N):
            xbar += xs[i]
        xbar /= N
            
        for i in range(N):
            cdif = xs[i] - xbar
            norm += cdif*cdif
            
            if i < N-k:
                rho += cdif * (xs[i+k]-xbar)
     
        rho *= N / ((N-k)*norm)
    
    return <_fptype> rho
        
cpdef autocorrel_elem(xs, k):
    '''Calculate the ``k``-lag sample autocorrelation of ``xs``. ``xs`` must be float32 or float64.'''
        
    cdef:
        int typecode = PyArray_TYPE(xs)
    
    if typecode == NPY_FLOAT64:
        return _autocorrel_elem[numpy.float64_t](xs,k)
    elif typecode == NPY_FLOAT32:
        return _autocorrel_elem[numpy.float32_t](xs,k)
    else:
        raise TypeError('unsupported type')
    
cpdef Py_ssize_t get_bssize(double alpha) nogil:
    '''Return a bootstrap data set size appropriate for the given confidence level.'''
    
    cdef:
        Py_ssize_t bssize=1, i
        
    #return int(10**(math.ceil(-math.log10(alpha)) + 1))
    for i in range(<long> ceil(-log10(alpha)) + 1):
        bssize *= 10
    return bssize

cpdef mcbs_ci(dataset, estimator, alpha, dlen, n_sets=None, args=None, kwargs=None, sort=numpy.msort):
    '''Perform a Monte Carlo bootstrap estimate for the (1-``alpha``) confidence interval
    on the given ``dataset`` with the given ``estimator``.  This routine is not appropriate
    for time-correlated data.
    
    Returns ``(estimate, ci_lb, ci_ub)`` where ``estimate`` is the application of the
    given ``estimator`` to the input ``dataset``, and ``ci_lb`` and ``ci_ub`` are the
    lower and upper limits, respectively, of the (1-``alpha``) confidence interval on 
    ``estimate``. 
    
    ``estimator`` is called as ``estimator(dataset, *args, **kwargs)``. Common estimators include:
      * numpy.mean -- calculate the confidence interval on the mean of ``dataset``
      * numpy.median -- calculate a confidence interval on the median of ``dataset``
      * numpy.std -- calculate a confidence interval on the standard deviation of ``datset``.

    ``n_sets`` is the number of synthetic data sets to generate using the given ``estimator``,
    which will be chosen using `get_bssize()`_ if ``n_sets`` is not given. 
    
    ``sort`` can be used
    to override the sorting routine used to calculate the confidence interval, which should
    only be necessary for estimators returning vectors rather than scalars.    
    '''
    
    if alpha > 0.5:
        raise ValueError('alpha ({}) > 0.5'.format(alpha))
    
    args = args or ()
    kwargs = kwargs or {}
    
    # dataset SHOULD be a dictionary.
    d_input = dataset.copy()
    # Here, we're dumping in any extra kwarg arguments to pass in to the estimator.
    try:
        d_input.update(kwargs)
    except:
        pass

    fhat = estimator(**d_input)
    
    try:
        estimator_shape = fhat.shape
    except AttributeError:
        estimator_shape = ()
        
    try:
        estimator_dtype = fhat.dtype
    except AttributeError:
        estimator_dtype = type(fhat)
        
    n_sets = n_sets or get_bssize(alpha)

    f_synth = numpy.empty((n_sets,) + estimator_shape, dtype=estimator_dtype)
    
    for i in xrange(n_sets):
        indices = numpy.random.randint(dlen, size=(dlen,))
        d_synth = {}
        for key, dset in dataset.iteritems():
            d_synth[key] = numpy.take(dset, indices, axis=0)
        d_input = d_synth.copy()
        try:
            d_input.update(kwargs)
        except:
            pass
        f_synth[i] = estimator(**d_input)
        del indices
        
    f_synth_sorted = sort(f_synth)
    lbi = int(math.floor(n_sets*alpha/2.0))
    ubi = int(math.ceil(n_sets*(1-alpha/2.0)))                     
    lb = f_synth_sorted[lbi]
    ub = f_synth_sorted[ubi]
    sterr = numpy.std(f_synth_sorted)
    
    del f_synth_sorted, f_synth
    return (fhat, lb, ub, sterr)

cpdef mcbs_correltime(dataset, alpha, n_sets = None):
    '''Calculate the correlation time of the given ``dataset``, significant to the
    (1-``alpha``) level, using the method described in Huber & Kim, "Weighted-ensemble
    Brownian dynamics simulations for protein association reactions" (1996), 
    doi:10.1016/S0006-3495(96)79552-8. An appropriate balance between space and speed
    is chosen based on the size of the input data.
    
    Returns 0 for data statistically uncorrelated with (1-alpha) confidence, otherwise
    the correlation length. (Thus, the appropriate stride for blocking is the 
    result of this function plus one.)'''
    
    n_sets = n_sets or get_bssize(alpha)
    if dataset.nbytes * n_sets > CORRELTIME_CROSSOVER:
        return mcbs_correltime_fullauto(dataset, alpha, n_sets)
    else:
        return mcbs_correltime_small(dataset, alpha, n_sets)

cpdef mcbs_correltime_fullauto(dataset, alpha, n_sets):
    '''Specialization of `mcbs_correltime`_ for large datasets, where rather than generating
    ``n_sets`` synthetic datasets and then evaluating autocorrelation elements on them, a
    synthetic data set is generated for each autocorrelation element required. This trades time
    for space.'''
    
    for k in xrange(1,len(dataset)//2):        
        data_rho, synrho_lb, synrho_ub = mcbs_ci(dataset, autocorrel_elem, alpha, n_sets, args=(k,))
        if data_rho > synrho_lb and data_rho < synrho_ub:
            return k-1
    else:
        return len(dataset)
        #raise ValueError('correlation length exceeds half the data set length')

cpdef mcbs_correltime_small(dataset, alpha, n_sets):
    '''Specialization of `mcbs_correltime`_ for small-to-mediam datasets, where ``n_sets``
    synthetic datasets are generated and stored, and then autocorrelation elements
    calculated on them. This implementation trades space for time.'''
    
    if alpha > 0.5:
        raise ValueError('alpha ({}) > 0.5'.format(alpha))    
    
    dlen = len(dataset)
    lbi = int(math.floor(n_sets*alpha/2.0))
    ubi = int(math.ceil(n_sets*(1-alpha/2.0)))
    
    synth_sets = numpy.empty((n_sets,)+dataset.shape, dtype=dataset.dtype)
    synth_acf_elems = numpy.empty((n_sets,), numpy.float64)
    
    for i in xrange(n_sets):
        synth_sets[i] = numpy.take(dataset,numpy.random.randint(dlen,size=(dlen,)))
    
    for k in xrange(1,dlen//2):
        data_rho = autocorrel_elem(dataset, k)
        
        for i in xrange(n_sets):
            synth_acf_elems[i] = autocorrel_elem(synth_sets[i],k)
            
        synth_acf_elems.sort()
        
        if data_rho > synth_acf_elems[lbi] and data_rho < synth_acf_elems[ubi]:
            return k-1
    else:
        #raise ValueError('correlation length exceeds half the data set length')
        return len(dataset)


