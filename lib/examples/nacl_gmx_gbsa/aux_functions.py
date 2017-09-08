#!/usr/bin/env python
#
# aux_functions.py
#
# This Python module defines auxilliary functions used during the WESTPA 
# simulation. In particular, these functions load atomic coordinate data
# and log files that are output by GROMACS.
#

import numpy

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

    # Parse data
    line_i  = 0
    frame_i = 0
    field_i = 0
    while line_i < len(raw_text):
        line = raw_text[line_i]
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
