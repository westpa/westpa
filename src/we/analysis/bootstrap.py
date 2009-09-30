from __future__ import division
import numpy
from math import floor, ceil

def mc_bootstrap(data, n_sets):
    """Construct n_sets of synthetic data by taking random rows from data.
    Returns a list of the synthetic data sets"""
    synth_sets = []
    N = data.shape[0]
    for i in xrange(0, n_sets):
        indices = numpy.random.randint(N, size=(N,))
        synth = numpy.empty_like(data)
        for r in xrange(0, N):
            synth[r] = data[indices[r]]
        synth_sets.append(synth)
    return synth_sets

def mc_bootstrap_a(data, n_sets):
    """Construct n_sets of synthetic data by taking random rows from data.
    Returns a new array whose first dimension indexes the synthetic data sets
    and whose subsequent indices match those of data"""
    synth = numpy.empty((n_sets,) + data.shape, data.dtype)
    N = data.shape[0]
    for i in xrange(0, n_sets):
        indices = numpy.random.randint(N, size=(N,))
        synth[i] = data[indices]
    return synth

def bootstrap_ci(estimator, data, alpha, n_sets=1000):
    """Construct a confidence interval for the given estimator of the given
    data using Monte Carlo bootstrap.  Returns (estimated_value, ci_lb, ci_ub) 
    where estimated_value = estimator(data) and ci_lb and ci_ub are the lower 
    and upper bounds of the (1-alpha)% confidence interval, respectively.
    """
    
    # Santize alpha, just in case someone uses 0.95 instead of 0.05
    alpha = min(alpha, abs(1-alpha))
    half_alpha = alpha / 2.0 
    
    fhat = estimator(data)
    f_synth = numpy.empty((n_sets,), data.dtype)
    synth = mc_bootstrap_a(data, n_sets)
    for i in xrange(0, n_sets):
        f_synth[i] = estimator(synth[i])
    f_synth.sort()
    lb = f_synth[int(floor(n_sets*half_alpha))]
    ub = f_synth[int(ceil(n_sets*(1-half_alpha)))]
    
    return (fhat, lb, ub)
