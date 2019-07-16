import os, sys, math, itertools
import numpy
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

class System(WESTSystem):
    """
    Over-damped Langevin dynamics
    """

    def initialize(self):
        self.pcoord_ndim = 1
        self.pcoord_len = 21
        self.pcoord_dtype = numpy.float32

        self.bin_mapper = RectilinearBinMapper([ list(numpy.arange(0.0, 10.1, 0.1)) ])
        self.bin_target_counts = numpy.empty((self.bin_mapper.nbins,), numpy.int_)
        self.bin_target_counts[...] = 10

def displacement_loader(fieldname, coord_filename, segment, single_point=False):
    """
    Loads and stores coordinates

    **Arguments:**
        :*fieldname*:      Key at which to store dataset
        :*coord_filename*: Temporary file from which to load coordinates
        :*segment*:        WEST segment
        :*single_point*:   Data to be stored for a single frame
                           (should always be false)
    """

    #print('kfw inside coord_loader', coord_filename)

    # Load coordinates
    n_frames = 21
    coord    = numpy.loadtxt(coord_filename, dtype = numpy.float32)
    coord    = numpy.reshape(coord, (n_frames))

    # Save to hdf5
    segment.data[fieldname] = coord

def seg_status_loader(fieldname, log_filename, segment, single_point=False):
    """
    Loads and stores log

    **Arguments:**
        :*fieldname*:    Key at which to store dataset
        :*log_filename*: Temporary file from which to load log
        :*segment*:      WEST segment
        :*single_point*: Data to be stored for a single frame
                         (should always be false)
    """

    # Load log
#   with open(log_filename, 'r') as log_file:
#       raw_text = [line.strip() for line in log_file.readlines()]

    # Save to hdf5

    f = open(log_filename, 'r')
    seg_status = f.read()
    f.close 

    with open(log_filename, 'r') as log_file:
        for line in log_file:
            segment.status = int(line)

    #print('kfw inside seg_status_loader :', repr(segment.status))
