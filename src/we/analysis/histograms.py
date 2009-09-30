from __future__ import division
from copy import copy
import math, numpy
import bootstrap, probability

def dd_bins(data, bin_scheme=None):
    """Return a set of bins appropriate for a histogram of the given data.
    If scheme is 'scott', use D.W. Scott's choice for bin width, based on
    the sample standard deviation of the data set (see Biometrika 66(3), 
    605-610 [1979]). If scheme is 'freedman-diaconis' (or just 'freedman') use
    the Freedman-Diaconis choice for bin width, based on the interquartile
    range of the data set (see Z. Wahrschenlichkeitstheorie verw. Gebiete 57, 
    453-476 [1981]).
    
    Returns a sequence of length k+1, where k is the number of bins in the
    histogram.
    """
    
    scheme = bin_scheme or 'scott'
    nscheme = scheme.lower()
    if nscheme == 'scott':
        h = 3.5 * data.std() / len(data)**(1/3)
        mind = min(data)
        maxd = max(data)
    elif nscheme in ('freedman-diaconis', 'freedman'):
        dsort = data[data.argsort()]
        iqr = (dsort[int(math.ceil(-0.25*len(data)))] 
               - dsort[int(math.floor(0.25*len(data)))])
        h = 2*iqr / len(data)**(1/3)
        mind = dsort[0]
        maxd = dsort[-1]
    else:
        raise ValueError('invalid binning scheme %r' % scheme)
    return [mind + h*i for i in xrange(0, int(math.ceil((maxd-mind)/h)))]

def histogram_with_mc_errorbars(data, bin_scheme = None, 
                                alpha = 0.05, n_synth = 1000):
    """This is probably a bad idea.  Use a Q-Q plot to compare distributions
    instead, or use a Kolmogorov-Smirnov confidence band on an EDF"""
    
    bins = dd_bins(data, bin_scheme)
    
    (hist, hbins) = numpy.histogram(data, bins, normed=True, new=True)
    hist_lb = numpy.empty_like(hist)
    hist_ub = numpy.empty_like(hist)
    
    synth_hists = numpy.empty((n_synth,) + hist.shape, hist.dtype)
    synth_data = bootstrap.mc_bootstrap_a(data, n_synth)
    for (i, sd) in enumerate(synth_data):
        synth_hists[i,:] = (numpy.histogram(sd, bins, normed=True, new=True))[0]
        
    for ibin in xrange(0, synth_hists.shape[1]):
        hist_lb[ibin], hist_ub[ibin] = probability.ci_range(synth_hists[:,ibin], alpha)
    
    return (hist, hist_lb, hist_ub, bins)
