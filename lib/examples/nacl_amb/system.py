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
        self.pcoord_len   = 5
        self.pcoord_dtype = numpy.float32
        binbounds         = [ 0.00, 2.80, 2.88, 3.00, 3.10, 3.29, 3.79, 3.94,
                              4.12, 4.39, 5.43, 5.90, 6.90, 7.90, 8.90, 9.90,
                             10.90,11.90,12.90,13.90,14.90,15.90,float('inf')]
        self.bin_mapper   = RectilinearBinMapper([binbounds])

        self.bin_target_counts      = numpy.empty((self.bin_mapper.nbins,),
                                        numpy.int)
        self.bin_target_counts[...] = 24

def coord_loader(fieldname, temp_filename, segment, single_point=False):
    """
    Loads and stores coordinates

    **Arguments:**
        :*fieldname*:     Key at which to store dataset
        :*temp_filename*: Temporary file from which to load filename
        :*segment*:       WEST segment
        :*single_point*:  Data to be stored for a single frame
                          (should always be false)
    """
    import netCDF4

    # Load coordinates
    with open(temp_filename, 'r') as temp_file:
        coord_filename = temp_file.readlines()[0].strip()
    coord_file = netCDF4.Dataset(coord_filename)
    coord      = numpy.array(coord_file.variables["coordinates"])

    # Save to hdf5
    segment.data[fieldname] = coord

def log_loader(fieldname, temp_filename, segment, single_point=False):
    """
    Loads and stores log

    **Arguments:**
        :*fieldname*:     Key at which to store dataset
        :*temp_filename*: Temporary file from which to load filename
        :*segment*:       WEST segment
        :*single_point*:  Data to be stored for a single frame
                          (should always be false)
    """
    # Load log
    with open(temp_filename, 'r') as temp_file:
        log_filename = temp_file.readlines()[0].strip()
    with open(log_filename, 'r') as log_file:
        raw_text = [line.strip() for line in log_file.readlines()]

    # Determine number of frames and fields
    n_frames = 0
    n_fields = 0
    i        = 0
    while i < len(raw_text):
        if   raw_text[i].startswith("A V E R A G E S"):
            break
        elif raw_text[i].startswith("NSTEP"):
            n_frames += 1
            if n_fields == 0:
                while True:
                    if raw_text[i].startswith("----------"):
                        break
                    n_fields += len(raw_text[i].split("="))
                    i        += 1
        i += 1
    dataset = numpy.zeros((n_frames, n_fields), numpy.float32)

    # Parse data
    i = 0
    a = 0
    while i < len(raw_text):
        if   raw_text[i].startswith("A V E R A G E S"):
            break
        elif raw_text[i].startswith("NSTEP"):
            b = 0
            while True:
                if raw_text[i].startswith("----------"):
                    break
                line = raw_text[i].split("=")
                for j, key in enumerate(line[:-1]):
                    dataset[a,b] = line[j+1].split()[0]
                    b += 1
                i += 1
            a += 1
        i += 1

    # Save to hdf5
    if segment.initpoint_type == 2:
        dataset = dataset[1:]
    segment.data[fieldname] = dataset
