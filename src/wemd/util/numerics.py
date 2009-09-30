from __future__ import division
import numpy

def histogram_with_errorbars(data, bins = None, n_synth_sets = 1000):
    dmin = data.min()
    dmax = data.max()
    dstd = data.std()
    
    binwidth = 3.5*dstd/len(data)**(1/3)
    nbins = math.ceil((dmax-dmin)/binwidth)
    
    (hist, bins) = numpy.histogram(data, nbins, normed=True, new=True)
    hist_sd = numpy.empty_like(hist)
    
    synth_hists = numpy.empty((hist.shape[0], n_synth_sets), hist.dtype)
    synth_data = mc_bootstrap(data, n_synth_sets)
    for (i, sd) in enumerate(synth_data):
        (synth_hists[:,i], jbins) = numpy.histogram(sd, bins, normed=True, new=True)
        
    for i in xrange(0, hist_sd.shape[0]):
        hist_sd[i] = synth_hists[i].std()
        
    return (hist, hist_sd, bins)
