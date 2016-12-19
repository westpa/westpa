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

from _fasthist import histnd #@UnresolvedImport

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
