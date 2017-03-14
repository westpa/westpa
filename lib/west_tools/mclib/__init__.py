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

'''A package for performing Monte Carlo bootstrap estimates of
statistics.'''

from __future__ import print_function, division
import math, numpy

from _mclib import autocorrel_elem, mcbs_correltime, get_bssize, mcbs_ci #@UnresolvedImport

def mcbs_ci_correl(estimator_datasets, estimator, alpha, n_sets=None, args=None,
                   autocorrel_alpha = None, autocorrel_n_sets=None, subsample=None, 
                   do_correl=True, mcbs_enable=None, estimator_kwargs={}):
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
    if type(estimator_datasets).__name__ != 'dict':
        # Enforcing the data structure.
        pre_calculated = estimator_datasets
        estimator_datasets = {'a' : estimator_datasets}
        # This also probably means our estimator isn't going to handle kwargs, so we'll watch out for that later in testing.
        # We may have to replace the 'simple' estimator with a slightly more complex lambda function which simply ditches extra arguments.
    for key, dset in estimator_datasets.iteritems():
        estimator_datasets[key] = numpy.asanyarray(dset)
        dlen = dset.shape[0]

    # Why do we have 'estimator_datasets'?
    # Estimators may require many different sets of data to properly function; while we can send this in via the kwargs,
    # we may wish to decimate only a certain subset (due to the block bootstrapping) of the input parameters.
    # Therefore, 'estimator_datasets' should consist of datasets that must be sliced/decimated with the subsampling function.
    # Some estimators (such as the reweighting) may not be able to be decimated in a straightforward manner with a subsample function,
    # as we cannot pre-estimate the quantity without introducing error or bias.  In those cases, we may wish to pass on all the data,
    # but ensure that our estimator only includes certain iterations (and only in a certain way).
    
    n_sets = n_sets or get_bssize(alpha)
    autocorrel_n_sets = autocorrel_n_sets or get_bssize(autocorrel_alpha)

    if mcbs_enable == False:
        # While it's odd to support NOT doing the bootstrap in a library specifically designed for bootstrapping, 
        # supporting this functionality here makes writing the code a lot easier, as we can just pass in a flag.
        # Specifically, this is for situations in which error is not desired (that is, only a reasonable mean is desired).
        # It's often useful when doing a quick analysis.
        estimator_datasets.update(estimator_kwargs)
        try:
            estimator_datasets.update( { 'stride': 1 } )
        except:
            pass
        
        return_set = estimator(**estimator_datasets)
        # We don't try and pretend we're doing any error analysis.
        return return_set, return_set, return_set, 0, 1

    # We need to pre-generate the data; why not do it here?  We're already set up for it...
    precalc_kwargs = estimator_kwargs.copy()
    precalc_kwargs['stride'] = 1
    pre_calculated = []
    for block in range(1, dlen+1):
        for key, dset in estimator_datasets.iteritems():
            precalc_kwargs[key] = dset[0:block]
        pre_calculated.append(estimator(**precalc_kwargs))
    # We need to get rid of any NaNs.
    pre_calculated = numpy.asanyarray(pre_calculated)
    pre_calculated = pre_calculated[numpy.isfinite(pre_calculated)]
    # If this happens, we have a huge NaN problem.  That is, our estimator is failing to return meaningful
    # numbers.  We should catch this when it happens, and so raise an exception, here.
    # This is almost certainly due to estimator failure.  Double check that calculation.
    if pre_calculated.shape == (0,):
        raise NameError("Looks like the estimator failed.  This is likely a programming issue, and should be reported.")
    # If pre-calculated is not None, we'll use that instead of dataset.
    # We can also assume that it's a 1 dimensional set with nothing needed, so 'key' should work.
    if do_correl == True:
        correl_len = mcbs_correltime(pre_calculated, autocorrel_alpha, autocorrel_n_sets)
    else:
        correl_len = 0
    if correl_len == len(pre_calculated):
        # too correlated for meaningful calculations
        estimator_datasets.update(estimator_kwargs)
        try:
            estimator_datasets.update( { 'stride': 1 } )
        except:
            pass

        return estimator(**estimator_datasets), pre_calculated.min(), pre_calculated.max(), (numpy.std(pre_calculated)), correl_len

    # else, do a blocked bootstrap
    stride = correl_len + 1
    
    if stride == 1:
        # Some estimators may require the stride, so we pass it in.
        estimator_kwargs['stride'] = stride
        return mcbs_ci(dataset=estimator_datasets, estimator=estimator, alpha=alpha, dlen=dlen, n_sets=n_sets, args=args, kwargs=estimator_kwargs, sort=numpy.msort) + (correl_len,)
    else:
        subsample = subsample or (lambda x: x[numpy.random.randint(len(x))])
        # Let's make sure we decimate every array properly...
        decim_list = {}
        for key,dset in estimator_datasets.iteritems():
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
            estimator_kwargs['stride'] = stride
        
        return mcbs_ci(dataset=decim_list, estimator=estimator, alpha=alpha, dlen=dlen, n_sets=n_sets, args=args, kwargs=estimator_kwargs, sort=numpy.msort) + (correl_len,)


# These are blocks designed to evaluate simple information sets.
# Whether they should go here or in westtoools is somewhat up for debate.
# Currently, nothing actually uses them, so there's that.

def _1D_simple_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, mcbs_enable, subsample=numpy.mean, **extra):
    # This is actually appropriate for anything with a directly measured, 1D dataset, i.e.,
    # Fluxes, color populations, and state populations.
    results = []
    for istate in xrange(nstates):
        # Not sure if we need a jstate for these estimators, but we'll see.
        kwargs = { 'istate' : istate , 'jstate': 'B'}
        estimator_datasets = {'dataset': data_input['dataset'][:,istate]}
        ci_res = mcbs_ci_correl(estimator_datasets,estimator=(lambda stride, dataset: numpy.mean(dataset)),
                                alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                subsample=subsample, do_correl=do_correl, mcbs_enable=mcbs_enable)

        results.append((name, iblock,istate,(start,stop)+ci_res))

    return results

def _2D_simple_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, mcbs_enable, subsample=numpy.mean, **extra):
    # This is really just a simple 2D block for less complex datasets, but there it is.
    # It's probably limited in this use case to conditional_fluxes, but anything that's an i to j process that is directly measured
    # is suitable for use with this.
    results = []
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            if istate == jstate: continue
            kwargs = { 'istate' : istate, 'jstate': jstate }
            #dataset = {'dataset': cond_fluxes[:, istate, jstate]}
            estimator_datasets = {'dataset': data_input['dataset'][:, istate, jstate] }
            ci_res = mcbs_ci_correl(estimator_datasets,estimator=(lambda stride, dataset: numpy.mean(dataset)),
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=subsample, do_correl=do_correl, mcbs_enable=mcbs_enable)

            results.append((name, iblock, istate, jstate, (start,stop) + ci_res))

    return results
