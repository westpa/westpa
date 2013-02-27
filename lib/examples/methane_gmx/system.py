from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

class System(WESTSystem):
    def initialize(self):
        self.pcoord_ndim = 1
        self.pcoord_len = 51
        self.pcoord_dtype = numpy.float32
        binbounds = [0.0] + [4.0+1.0*i for i in xrange(0,13)] + [30,float('inf')]
        self.bin_mapper = RectilinearBinMapper([binbounds])
        self.bin_target_counts = numpy.empty((self.bin_mapper.nbins,), numpy.int)
        self.bin_target_counts[...] = 48

def coord_loader(fieldname, coord_file, segment, single_point=False):
    coord_raw = numpy.loadtxt(coord_file, dtype=numpy.float32) 

    npts = len(coord_raw)
    assert coord_raw.shape[1] % 3 == 0
    ngrps = coord_raw.shape[1] // 3

    coords = numpy.empty((ngrps, npts, 3), numpy.float32)
    for igroup in xrange(ngrps):
        for idim in xrange(3):
            coords[igroup,:,idim] = coord_raw[:,igroup*3+idim]
    # convert to Angstroms
    coords *= 10

    segment.data[fieldname] = coords
    

