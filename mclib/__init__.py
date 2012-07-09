from __future__ import print_function, division
from _mclib import autocorrel_elem

import math, numpy

def get_bssize(alpha):
    '''Return a bootstrap data set size appropriate for the given confidence level'''
    return int(10**(math.ceil(-math.log10(alpha)) + 1))

def mcbs_ci(dataset, estimator, alpha, n_sets=None, args=None, kwargs=None, sort=numpy.msort):
    
    if alpha > 0.5:
        raise ValueError('alpha ({}) > 0.5'.format(alpha))
    
    args = args or ()
    kwargs = kwargs or {}
    dataset = numpy.asanyarray(dataset)
    dlen = len(dataset)
    
    fhat = estimator(dataset, *args, **kwargs)
    
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
        f_synth[i] = estimator(numpy.take(dataset,indices), *args, **kwargs)
        del indices
        
    f_synth_sorted = sort(f_synth)
    lbi = int(math.floor(n_sets*alpha/2.0))
    ubi = int(math.ceil(n_sets*(1-alpha/2.0)))                     
    lb = f_synth_sorted[lbi]
    ub = f_synth_sorted[ubi]
    
    del f_synth_sorted, f_synth
    return (fhat, lb, ub)

def mcbs_correltime(dataset, alpha, n_sets=None):
    dlen = len(dataset)
    n_sets = n_sets or get_bssize(alpha)
    
    for k in xrange(1,len(dataset)//2):        
        data_rho, synrho_lb, synrho_ub = mcbs_ci(dataset, autocorrel_elem, alpha, n_sets, args=(k,))
        if data_rho > synrho_lb and data_rho < synrho_ub:
            return k-1
    else:
        raise ValueError('correlation length exceeds half the data set length')

def mcbs_correltime_small(dataset, alpha, n_sets=None):
    dlen = len(dataset)
    n_sets = n_sets or get_bssize(alpha)
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
        raise ValueError('correlation length exceeds half the data set length')

def mcbs_ci_correl(dataset, estimator, alpha, n_sets=None, args=None, kwargs=None, sort=numpy.msort,
                   decorrel=None):
    
    if alpha > 0.5:
        raise ValueError('alpha ({}) > 0.5'.format(alpha))
    
    if dataset.ndim != 1:
        raise ValueError('correlated time series MCBS analysis available only for 1-dimensional data')
    
    dataset = numpy.asanyarray(dataset)
    dlen = len(dataset)
    n_sets = n_sets or get_bssize(alpha)
    
    correl_len = mcbs_correltime_small(dataset, alpha, n_sets)
    stride = correl_len + 1
    
    if stride == 1:
        return mcbs_ci(dataset, estimator, alpha, n_sets, args, kwargs, sort)
    else:
        n_slices = dlen // stride
        #if dlen % stride > 0: n_slices += 1
        decorrel = decorrel or (lambda x: x[numpy.random.randint(len(x))])
        decorrel_set = numpy.empty((n_slices,), dtype=dataset.dtype)
        for iout, istart in enumerate(xrange(0,dlen-stride+1,stride)):
            sl = dataset[istart:istart+stride]
            decorrel_set[iout] = decorrel(sl)
        
        return mcbs_ci(decorrel_set, estimator, alpha, n_sets, args, kwargs, sort)
