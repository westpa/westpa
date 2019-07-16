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


'''
Tools for Monte Carlo bootstrap error analysis
'''


import logging

log = logging.getLogger(__name__)

import math, numpy

import westpa
from oldtools.aframe import AnalysisMixin

class MCBSMixin(AnalysisMixin):
    def __init__(self):
        super(MCBSMixin,self).__init__()
        self.mcbs_alpha = None
        self.mcbs_nsets = None
        self.mcbs_display_confidence = None
        
    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc=super(MCBSMixin,self).add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
        group = parser.add_argument_group('Monte Carlo bootstrap options')
        group.add_argument('--confidence', dest='mcbs_confidence', type=float, default=0.95, metavar='P',
                           help='''Construct a confidence interval of width P (default: 0.95=95%%).''')
        group.add_argument('--bssize', dest='mcbs_nsets', type=int, metavar='NSETS',
                           help='''Use NSETS synthetic data sets to calculate confidence intervals (default:
                                   calculated based on confidence level, but not less than 1000).''' )
        
    def process_args(self, args, upcall = True):
        self.mcbs_alpha = 1-args.mcbs_confidence
        self.mcbs_nsets = args.mcbs_size if args.mcbs_nsets else min(1000,calc_mcbs_nsets(self.mcbs_alpha))
        self.mcbs_display_confidence = '{:.{cp}f}'.format(100*args.mcbs_confidence,
                                                                    cp = -int(math.floor(math.log10(self.mcbs_alpha)))-2) 
        westpa.rc.pstatus('Using bootstrap of {:d} sets to calculate {:s}% confidence interval (alpha={:g}).'
                        .format(self.mcbs_nsets, self.mcbs_display_confidence, self.mcbs_alpha))

        if upcall:
            try:
                upfunc = super(MCBSMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)
                
    def calc_mcbs_nsets(self, alpha = None):
        alpha = alpha or self.mcbs_alpha
        return calc_mcbs_nsets(alpha)

    def calc_ci_bound_indices(self, n_sets=None, alpha=None):
        n_sets = n_sets or self.mcbs_nsets
        alpha = alpha or self.mcbs_alpha
        return calc_ci_bound_indices(n_sets, alpha)
            

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])
                
def calc_mcbs_nsets(alpha):
    '''Return a bootstrap data set size appropriate for the given confidence level.'''
    return int(10**(math.ceil(-math.log10(alpha))+1))

def calc_ci_bound_indices(n_sets, alpha):
    return (int(math.floor(n_sets*alpha/2)), int(math.ceil(n_sets*(1-alpha/2))))

def bootstrap_ci_ll(estimator, data, alpha, n_sets, storage, sort, eargs=(), ekwargs={}, fhat=None):
    '''Low-level routine for calculating bootstrap error estimates.  Arguments and return values are as those for
    ``bootstrap_ci``, except that no argument is optional except additional arguments for the estimator (``eargs``, ``ekwargs``).
    ``data`` must be an array (or subclass), and an additional array ``storage`` must be provided, which
    must be appropriately shaped and typed to hold ``n_sets`` results from ``estimator``.  Further, if the
    value ``fhat`` of the estimator must be pre-calculated to allocate ``storage``, then its value may be
    passed; otherwise, ``estimator(data,*eargs,**kwargs)`` will be called to calculate it.'''
    
    if fhat is None:
        fhat = estimator(data, *eargs, **ekwargs)
    dlen = len(data)
    
    for iset in range(n_sets):
        indices = numpy.random.randint(dlen,size=(dlen,))
        storage[iset] = estimator(data[indices], *eargs, **ekwargs)
        
    synth_sorted = sort(storage)
    lbi = int(math.floor(n_sets*alpha/2))
    ubi = int(math.ceil(n_sets*(1-alpha/2)))
    
    lb = synth_sorted[lbi]
    ub = synth_sorted[ubi]
    
    try:
        return (fhat, lb, ub, ub-lb, abs((ub-lb)/fhat) if fhat else 0, max(ub-fhat,fhat-lb))
    finally:
        del fhat, lb, ub, indices

def bootstrap_ci(estimator, data, alpha, n_sets=None, sort=numpy.msort, eargs=(), ekwargs={}):
    '''Perform a Monte Carlo bootstrap of a (1-alpha) confidence interval for the given ``estimator``.
    Returns (fhat, ci_lower, ci_upper), where fhat is the result of ``estimator(data, *eargs, **ekwargs)``,
    and ``ci_lower`` and ``ci_upper`` are the lower and upper bounds of the surrounding confidence
    interval, calculated by calling ``estimator(syndata, *eargs, **ekwargs)`` on each synthetic data
    set ``syndata``.  If ``n_sets`` is provided, that is the number of synthetic data sets generated,
    otherwise an appropriate size is selected automatically (see ``calc_mcbs_nsets()``).
    
    ``sort``, if given, is applied to sort the results of calling ``estimator`` on each 
    synthetic data set prior to obtaining the confidence interval. This function must sort
    on the last index.
    
    Individual entries in synthetic data sets are selected by the first index of ``data``, allowing this
    function to be used on arrays of multidimensional data.
    
    Returns (fhat, lb, ub, ub-lb, abs((ub-lb)/fhat), and max(ub-fhat,fhat-lb)) (that is, the estimated value, the
    lower and upper bounds of the confidence interval, the width of the confidence interval, the relative
    width of the confidence interval, and the symmetrized error bar of the confidence interval).'''

    data = numpy.asanyarray(data)
    fhat = numpy.squeeze(estimator(data, *eargs, **ekwargs))
    n_sets = n_sets or calc_mcbs_nsets(alpha)
    fsynth = numpy.empty((n_sets,), dtype=fhat.dtype)
    try:
        return bootstrap_ci_ll(estimator, data, alpha, n_sets or calc_mcbs_nsets(alpha), fsynth, sort, eargs, ekwargs, fhat)
    finally:
        del fsynth
