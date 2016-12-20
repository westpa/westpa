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

'''A package for performing Monte Carlo bootstrap estimates of
statistics.'''

from __future__ import print_function, division
import math, numpy

from _mclib import autocorrel_elem, mcbs_correltime, get_bssize, mcbs_ci #@UnresolvedImport

def mcbs_ci_correl(dataset, estimator, alpha, n_sets=None, args=None, kwargs=None,
                   autocorrel_alpha = None, autocorrel_n_sets=None, subsample=None, pops=None,
                   istate=None, jstate=None, correl=True):
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

    if kwargs == None:
        kwargs = {}
    # We're adding in this stuff.  Bit hackish, but.
    
    if correl == True:
        correl_len = mcbs_correltime(dataset, autocorrel_alpha, autocorrel_n_sets)
    else:
        correl_len = 0
        del kwargs['correl']
    if pops != None:
        kwargs['pops'] = pops
        kwargs['istate'] = istate
        kwargs['jstate'] = jstate
    if correl_len == len(dataset):
        #if pops == None:
        # too correlated for meaningful calculations
        return estimator(dataset, *(args or ()), **(kwargs or {})), dataset.min(), dataset.max(), correl_len
        #else:
        #    return estimator(dataset, dataset_two=pops, istate=istate, jstate=jstate, *(args or ()), **(kwargs or {})), dataset.min(), dataset.max(), correl_len
        
    # else, do a blocked bootstrap
    stride = correl_len + 1
    
    if stride == 1:
        return mcbs_ci(dataset, estimator, alpha, n_sets, args, kwargs, numpy.msort) + (correl_len,)
        #else:
            #return mcbs_ci(dataset, estimator, alpha, n_sets, args, kwargs, numpy.msort, dataset_two=pops, istate=istate, jstate=jstate) + (correl_len,)
    else:
        n_slices = dlen // stride
        #if dlen % stride > 0: n_slices += 1
        subsample = subsample or (lambda x: x[numpy.random.randint(len(x))])
        decim_set = numpy.empty((n_slices,), dtype=dataset.dtype)
        if pops != None:
            decim_pops = numpy.empty((n_slices, pops.shape[1]), dtype=dataset.dtype)
        for iout, istart in enumerate(xrange(0,dlen-stride+1,stride)):
            #decim_set[iout] = subsample(sl, **kwargs)
            if pops == None:
                sl = dataset[istart:istart+stride]
                decim_set[iout] = subsample(sl)
            else:
                sl = dataset[istart:istart+stride]
                decim_set[iout] = subsample(sl)
                pl = pops[istart:istart+stride, :]
                decim_pops[iout, :] = subsample(pl, axis=0)
                #decim_set[iout] /= n_slices
                #decim_pops[iout, :] /= n_slices
        if pops != None:
            kwargs['pops'] = decim_pops
        
        #return mcbs_ci(decim_set, estimator, alpha, n_sets, args, kwargs, numpy.msort, dataset_two=pops, istate=istate, jstate=jstate) + (correl_len,)
        return mcbs_ci(decim_set, estimator, alpha, n_sets, args, kwargs, numpy.msort) + (correl_len,)

def mcbs_ci_correl_rw(dataset, estimator, alpha, n_sets=None, args=None,
        autocorrel_alpha = None, autocorrel_n_sets=None, subsample=None, pre_calculated=None, correl=False, **kwargs):
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
    
    # We're now passing in dataset as a dict, so we need to enforce that for compatibility with older tools.
    # This just takes our dataset and puts it into a dict, as it's likely that we're using
    # mean or median as our estimators, which take "a" as argument input.
    if type(dataset).__name__ != 'dict':
        # Enforcing the data structure.
        pre_calculated = dataset
        dataset = {'a' : dataset}
    for key, dset in dataset.iteritems():
        dataset[key] = numpy.asanyarray(dset)
        dlen = dset.shape[0]
    
    n_sets = n_sets or get_bssize(alpha)
    autocorrel_n_sets = autocorrel_n_sets or get_bssize(autocorrel_alpha)

    # We probably need to get this from the rates, so we'll still end up passing those in.
    # If pre-calculated is not None, we'll use that instead of dataset.
    # We can also assume that it's a 1 dimensional set with nothing needed, so 'key' should work.
    if correl == True:
        correl_len = 0
    else:
        correl_len = mcbs_correltime(pre_calculated, autocorrel_alpha, autocorrel_n_sets)
    if correl_len == len(pre_calculated):
        # too correlated for meaningful calculations
        d_input = dataset.copy()
        kwargs['stride'] = 1
        try:
            d_input.update(kwargs)
        except:
            pass

        return estimator(**d_input), pre_calculated.min(), pre_calculated.max(), (numpy.std(pre_calculated)), correl_len
        
    # else, do a blocked bootstrap
    stride = correl_len + 1
    
    if stride == 1:
        #return mcbs_ci(dataset, estimator, alpha, dlen, n_sets, args, kwargs, numpy.msort) + (correl_len,)
        kwargs['stride'] = stride
        return mcbs_ci(dataset=dataset, estimator=estimator, alpha=alpha, dlen=dlen, n_sets=n_sets, args=args, kwargs=kwargs, sort=numpy.msort) + (correl_len,)
    else:
        #n_slices = dlen // stride
        #if dlen % stride > 0: n_slices += 1
        subsample = subsample or (lambda x: x[numpy.random.randint(len(x))])
        #decim_set = numpy.empty((n_slices,), dtype=dataset.dtype)
        # Let's make sure we decimate every array properly...
        decim_list = {}
        for key,dset in dataset.iteritems():
            dset_shape = list(dset.shape)
            n_slices = dset_shape[0] // stride
            dset_shape[0] = n_slices
            decim_set = numpy.empty((dset_shape), dtype=dset.dtype)
            for iout, istart in enumerate(xrange(0,dset.shape[0]-stride+1,stride)):
                sl = dset[istart:istart+stride]
                # We assume time is the 0th axis.
                # Okay, so non-optimal.  Population requires the axis subsampling to be done just so...
                try:
                    decim_set[iout] = subsample(sl, axis=0)
                except:
                    decim_set[iout] = subsample(sl)
            decim_list[key] = decim_set
            dlen = dset_shape[0]
            kwargs['stride'] = stride
        
        return mcbs_ci(dataset=decim_list, estimator=estimator, alpha=alpha, dlen=dlen, n_sets=n_sets, args=args, kwargs=kwargs, sort=numpy.msort) + (correl_len,)
