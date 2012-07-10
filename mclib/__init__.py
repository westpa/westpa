'''A package for performing Monte Carlo bootstrap estimates of
statistics.'''

from __future__ import print_function, division
import math, numpy

#class OvercorrelationError(ValueError):
#    '''Raised by bootstrapping routines for time-correlated data when
#    the correlation length is too large to take correlation into account
#    meaningfully.'''
#    pass

# Number of bytes in a dataset above which correlation time estimates are 
# not performed by generating synthetic data sets once and then evaluating 
# the autocorrelation elements on them, but rather by generating a new
# synthetic dataset each time an autocorrelation element is required.
# Essentially, this limits the RAM used during the correlation time
# calculation to approximately CORRELTIME_CROSSOVER.  Keep in mind that
# the space savings is at the cost of needing *many* more random numbers,
# which will almost certainly become the rate-limiting step of the
# calculation.
CORRELTIME_CROSSOVER = 512*1024*1024

try:
    # by default, try to use the accelerated autocorrelation element 
    # routines for 
    from _mclib import autocorrel_elem
except ImportError:
    def autocorrel_elem(xs, k):
        '''Calculate the ``k``-lag sample autocorrelation of ``xs``'''
        xsmm = xs - xs.mean()
        N = len(xs)
        return N/(N-k) * (xsmm[k:]*xsmm[:-k]).sum() / (xsmm*xsmm).sum()

def get_bssize(alpha):
    '''Return a bootstrap data set size appropriate for the given confidence level.'''
    return int(10**(math.ceil(-math.log10(alpha)) + 1))

def mcbs_ci(dataset, estimator, alpha, n_sets=None, args=None, kwargs=None, sort=numpy.msort):
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

def mcbs_correltime(dataset, alpha, n_sets = None):
    '''Calculate the correlation time of the given ``dataset``, significant to the
    (1-``alpha``) level, using the method described in Huber & Kim, "Weighted-ensemble
    Brownian dynamics simulations for protein association reactions" (1996), 
    doi:10.1016/S0006-3495(96)79552-8. An appropriate balance between space and speed
    is chosen based on the size of the input data.'''
    
    n_sets = n_sets or get_bssize(alpha)
    if dataset.nbytes * n_sets > CORRELTIME_CROSSOVER:
        return mcbs_correltime_fullauto(dataset, alpha, n_sets)
    else:
        return mcbs_correltime_small(dataset, alpha, n_sets)

def mcbs_correltime_fullauto(dataset, alpha, n_sets=None):
    '''Specialization of `mcbs_correltime`_ for large datasets, where rather than generating
    ``n_sets`` synthetic datasets and then evaluating autocorrelation elements on them, a
    synthetic data set is generated for each autocorrelation element required. This trades time
    for space.'''
    
    n_sets = n_sets or get_bssize(alpha)
    
    for k in xrange(1,len(dataset)//2):        
        data_rho, synrho_lb, synrho_ub = mcbs_ci(dataset, autocorrel_elem, alpha, n_sets, args=(k,))
        if data_rho > synrho_lb and data_rho < synrho_ub:
            return k-1
    else:
        return len(dataset)
        #raise ValueError('correlation length exceeds half the data set length')

def mcbs_correltime_small(dataset, alpha, n_sets=None):
    '''Specialization of `mcbs_correltime`_ for small-to-mediam datasets, where ``n_sets``
    synthetic datasets are generated and stored, and then autocorrelation elements
    calculated on them. This implementation trades space for time.'''
    
    if alpha > 0.5:
        raise ValueError('alpha ({}) > 0.5'.format(alpha))    
    
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
        #raise ValueError('correlation length exceeds half the data set length')
        return len(dataset)

def mcbs_ci_correl(dataset, estimator, alpha, n_sets=None, args=None, kwargs=None,
                   autocorrel_alpha = None, subsample=None):
    '''Perform a Monte Carlo bootstrap estimate for the (1-``alpha``) confidence interval
    on the given ``dataset`` with the given ``estimator``.  This routine is appropriate
    for time-correlated data, using the method described in Huber & Kim, "Weighted-ensemble
    Brownian dynamics simulations for protein association reactions" (1996), 
    doi:10.1016/S0006-3495(96)79552-8 to determine a statistically-significant correlation time
    and then reducing the dataset by a factor of that correlation time before running a "classic"
    Monte Carlo bootstrap.

    Returns ``(estimate, ci_lb, ci_ub, correl_time)`` where ``estimate`` is the application of the
    given ``estimator`` to the input ``dataset``, ``ci_lb`` and ``ci_ub`` are the
    lower and upper limits, respectively, of the (1-``alpha``) confidence interval on 
    ``estimate``, and ``correl_time`` is the correlation time of the dataset, significant to
    (1-``autocorrel_alpha``).
    
    ``estimator`` is called as ``estimator(dataset, *args, **kwargs)``. Common estimators include:
      * numpy.mean -- calculate the confidence interval on the mean of ``dataset``
      * numpy.median -- calculate a confidence interval on the median of ``dataset``
      * numpy.std -- calculate a confidence interval on the standard deviation of ``datset``.

    ``n_sets`` is the number of synthetic data sets to generate using the given ``estimator``,
    which will be chosen using `get_bssize()`_ if ``n_sets`` is not given.
        
    ``autocorrel_alpha`` (which defaults to ``alpha``) can be used to adjust the significance
    level of the autocorrelation calculation. Note that too high a significance level (too low an
    alpha) for evaluating the significance of autocorrelation values can result in a failure to
    detect correlation if the autocorrelation function is noisy.    
    
    The given ``subsample`` function is used, if provided, to subsample the dataset prior to running
    the full Monte Carlo bootstrap. If none is provided, then a random entry from each correlated
    block is used as the value for that block.  Other reasonable choices include ``numpy.mean``,
    ``numpy.median``, ``(lambda x: x[0])`` or ``(lambda x: x[-1])``.  In particular, using
    ``subsample=numpy.mean`` will converge to the block averaged mean and standard error,
    while accounting for any non-normality in the distribution of the mean.
    '''
    
    if alpha > 0.5:
        raise ValueError('alpha ({}) > 0.5'.format(alpha))
    
    autocorrel_alpha = alpha if not autocorrel_alpha else autocorrel_alpha
    
    if dataset.ndim != 1:
        raise ValueError('correlated time series MCBS analysis available only for 1-dimensional data')
    
    dataset = numpy.asanyarray(dataset)
    dlen = len(dataset)
    n_sets = n_sets or get_bssize(alpha)
    
    correl_len = mcbs_correltime(dataset, autocorrel_alpha, n_sets)
    if correl_len == len(dataset):
        # too correlated for meaningful calculations
        return estimator(dataset, *(args or ()), **(kwargs or {})), dataset.min(), dataset.max(), correl_len
        
    # else, do a blocked bootstrap
    stride = correl_len + 1
    
    if stride == 1:
        return mcbs_ci(dataset, estimator, alpha, n_sets, args, kwargs, numpy.msort) + (correl_len,)
    else:
        n_slices = dlen // stride
        #if dlen % stride > 0: n_slices += 1
        subsample = subsample or (lambda x: x[numpy.random.randint(len(x))])
        decim_set = numpy.empty((n_slices,), dtype=dataset.dtype)
        for iout, istart in enumerate(xrange(0,dlen-stride+1,stride)):
            sl = dataset[istart:istart+stride]
            decim_set[iout] = subsample(sl)
        
        return mcbs_ci(decim_set, estimator, alpha, n_sets, args, kwargs, numpy.msort) + (correl_len,)
