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
    """
    System for sodium chloride association
    """

    def initialize(self):
        """
        Initializes system
        """
        self.pcoord_ndim  = 1
        self.pcoord_len   = 6
        self.pcoord_dtype = numpy.float32
        binbounds         = [ 0.00, 2.80, 2.88, 3.00, 3.10, 3.29, 3.79, 3.94,
                              4.12, 4.39, 5.43, 5.90, 6.90, 7.90, 8.90, 9.90,
                             10.90,11.90,12.90,13.90,14.90,15.90,float('inf')]
        self.bin_mapper   = RectilinearBinMapper([binbounds])

        self.bin_target_counts      = numpy.empty((self.bin_mapper.nbins,),
                                        numpy.int)
        self.bin_target_counts[...] = 24

def coord_loader(fieldname, coord_filename, segment, single_point=False):
    """
    Loads and stores coordinates

    **Arguments:**
        :*fieldname*:      Key at which to store dataset
        :*coord_filename*: Temporary file from which to load coordinates
        :*segment*:        WEST segment
        :*single_point*:   Data to be stored for a single frame
                           (should always be false)
    """
    # Load coordinates
    n_frames = 6
    n_atoms  = 2
    coord    = numpy.loadtxt(coord_filename, dtype = numpy.float32)
    coord    = numpy.reshape(coord, (n_frames, n_atoms, 3))

    # Save to hdf5
    segment.data[fieldname] = coord

def log_loader(fieldname, log_filename, segment, single_point=False):
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
    with open(log_filename, 'r') as log_file:
        raw_text = [line.strip() for line in log_file.readlines()]

    # Determine number of fields
    n_frames = 6
    n_fields = 0
    line_i   = 0
    starts   = []
    while line_i < len(raw_text):
        line = raw_text[line_i]
        if len(line.split()) > 0:
            start = line.split()[0]
            if start in starts:
                break
            else:
                try:
                    float(start)
                    n_fields += len(line.split())
                except ValueError:
                    starts.append(start)
        line_i += 1
    dataset = numpy.zeros((n_frames, n_fields), numpy.float32)
#    print(dataset.shape, starts)

    # Parse data
    line_i  = 0
    frame_i = 0
    field_i = 0
    while line_i < len(raw_text):
        line = raw_text[line_i]
#        print(line_i, frame_i, field_i, line)
        if len(line.split()) > 0:
            start = line.split()[0]
            try:
                float(start)
            except ValueError:
                starts.append(start)
            if start not in starts:
                for field in line.split():
                    dataset[frame_i, field_i] = float(field)
                    if field_i == n_fields - 1:
                        frame_i += 1
                        field_i  = 0
                    else:
                        field_i += 1
        line_i += 1

    # Save to hdf5
    segment.data[fieldname] = dataset
