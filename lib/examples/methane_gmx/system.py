from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy
import west
from west import WESTSystem
from west.pcoords import PiecewiseRegionSet, ParticleSet, RectilinearRegionSet

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

class System(WESTSystem):

    def new_region_set(self):
        binbounds = [0.0] + [4.0+1.0*i for i in xrange(0,13)] + [30,float('inf')]
        region_set = RectilinearRegionSet([binbounds])
        for bin in region_set.get_all_bins():
            bin.target_count = 48
        return region_set

    def initialize(self):
        self.pcoord_ndim = 1
        self.pcoord_len = 51
        self.pcoord_dtype = numpy.float32

def coord_loader(fieldname, coord_file, segment, single_point=False):
    '''Load coordinates as output by g_traj and paste'''
    coord_raw = numpy.loadtxt(coord_file, dtype=numpy.float32)

    # each atom/group/C.O.M selected gets 4 columns:
    # time, x, y, z

    npts = len(coord_raw)
    assert coord_raw.shape[1] % 4 == 0
    ngrps = coord_raw.shape[1] // 4

    coords = numpy.empty((ngrps, npts, 3), numpy.float32)
    for igroup in xrange(ngrps):
        for idim in xrange(3):
            coords[igroup,:,idim] = coord_raw[:,igroup*4+idim+1]
    # convert to Angstroms
    coords *= 10

    segment.data[fieldname] = coords
    

