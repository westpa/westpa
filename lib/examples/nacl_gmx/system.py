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
    coord = numpy.loadtxt(coord_filename, dtype = numpy.float32)
    coord = numpy.reshape(coord, (int(coord.shape[0] / 2),
       int(coord.shape[1] / 3 * 2), 3))

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
    i        = 0
    n_frames = 0
    n_fields = 0
    while i < len(raw_text):
        if "A V E R A G E S" in raw_text[i]:
            break
        elif raw_text[i].startswith("Step"):
            n_frames += 1
            i        += 1
            if n_fields == 0:
                while True:
                    if raw_text[i].startswith("Step"):
                        i -= 1
                        break
                    elif (raw_text[i].startswith("Note")
                    or    raw_text[i].startswith("Grid")
                    or    raw_text[i].startswith("Writing")
                    or    raw_text[i].startswith("Large")):
                        pass
                    else:
                        for field in raw_text[i].split():
                            try:
                                float(field)
                                n_fields += 1
                            except:
                                pass
                    i += 1
        i += 1
    dataset = numpy.zeros((n_frames, n_fields), numpy.float32)

    # Parse data
    line_i  = 0
    frame_i = 0
    while line_i < len(raw_text):
        line = raw_text[line_i].split()
        if len(line) >= 1 and line[0] == "Step":
            field_i  = 0
            line_i  += 1
            while True:
                line = raw_text[line_i].split()
                if len(line) >= 1:
                    if line[0] == "Step":
                        frame_i += 1
                        field_i  = 0
                    elif line[0] == "<======":
                        line_i = len(raw_text)
                        break
                    else:
                        try:
                            float(line[0])
                            for field in line:
                                dataset[frame_i, field_i] = float(field)
                                field_i += 1
                        except ValueError:
                            pass
                line_i += 1
        line_i += 1

    # Save to hdf5
    segment.data[fieldname] = dataset
