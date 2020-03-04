
from ._fasthist import histnd

import numpy

def normhistnd(hist, binbounds):
    '''Normalize the N-dimensional histogram ``hist`` with corresponding
    bin boundaries ``binbounds``.  Modifies ``hist`` in place and returns
    the normalization factor used.'''

    ndim = hist.ndim
    
    if ndim != len(binbounds):
        raise ValueError('shape of histogram [{!r}] does not match bin boundary sets (there are {})'
                         .format(hist.shape, len(binbounds)))
    
    diffs = [numpy.diff(bb) for bb in binbounds]

    if ndim == 1:
        assert diffs[0].shape == hist.shape
        normfac = (hist * diffs[0]).sum()
    else:
        outers = numpy.multiply.outer(diffs[0], diffs[1])
        for delta in diffs[2:]:
            outers = numpy.multiply.outer(outers, delta)
        assert outers.shape == hist.shape, 'hist shape {} != outers shape {}'.format(hist.shape,outers.shape)
        # Divide by bin volumes
        hist /= outers
        normfac = hist.sum()
        #normfac = (hist * outers).sum()

    hist /= normfac
    return normfac
