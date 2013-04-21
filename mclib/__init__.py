# Copyright (C) 2013 Matthew C. Zwier
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

'''A package for performing Monte Carlo bootstrap estimates of
statistics.'''

from __future__ import print_function, division
import math, numpy

from _mclib import autocorrel_elem, mcbs_correltime, get_bssize, mcbs_ci #@UnresolvedImport

def mcbs_ci_correl(dataset, estimator, alpha, n_sets=None, args=None, kwargs=None,
                   autocorrel_alpha = None, autocorrel_n_sets=None, subsample=None):
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
    autocorrel_n_sets = autocorrel_n_sets or get_bssize(autocorrel_alpha)
    
    correl_len = mcbs_correltime(dataset, autocorrel_alpha, autocorrel_n_sets)
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
