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
        self.pcoord_len = 11
        self.pcoord_dtype = numpy.float32
        binbounds = ([0.0] + [2.8, 2.88, 3.0, 3.10, 3.29, 3.79, 3.94, 4.12, 4.39, 5.43] + [5.90+1.0*i for i in xrange(0,11)] + [30,float('inf')]) 
        self.bin_mapper = RectilinearBinMapper([binbounds])
        self.bin_target_counts = numpy.empty((self.bin_mapper.nbins,), numpy.int)
        self.bin_target_counts[...] = 24

def coord_loader(fieldname, coord_file, segment, single_point=False):
    coord_raw_file = open(coord_file, 'r') 
    coord_raw = coord_raw_file.readlines()
    coords = []
    for line in coord_raw:
        line = list([line[0:8], line[8:16], line[16:24], \
                     line[24:32], line[32:40], line[40:48]])
        line = map(numpy.float, line)
        line = numpy.array(line)
        coords.append(line) 
    coords = numpy.array(coords)
    assert coords.shape[1] % 3 == 0

    segment.data[fieldname] = coords
