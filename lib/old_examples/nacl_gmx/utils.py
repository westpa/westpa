#!/usr/bin/env python
import numpy

class TrajWriter(object):
    '''
    A class for writing out trajectory traces as an xyz file, for subsequent
    visualization.

    ---------
    Arguments
    ---------
    trace: A trace returned by w.trace.
    w: The w object loaded by w_ipa

    -----------------
    Keyword Arguments
    -----------------
    filename (default='trace.xyz'): A path designating the file to which the
      trajectory should be written, in xyz format.
    '''
    def __init__(self, trace, w, filename='trace.xyz'):
        self.trace = trace
        self.w = w
        self.filename = filename
        self._write()

    def _get_coords(self, iteration, seg_id):
        self.w.iteration = iteration
        coords = self.w.current.auxdata['coord'][seg_id]
        return coords

    def _write(self):
        all_coords = []
        starting_iteration = self.w.iteration 
        for i, iteration in enumerate(self.trace.iteration):
            seg_id = self.trace.seg_id[i]
            coords = self._get_coords(iteration, seg_id)
            # The last timepoint of one iteration is the same as the first
            # timepoint of the last, so skip the last timepoint of each 
            # iteration
            coords = coords[:-1]
            all_coords.append(coords)
        self.w.iteration = starting_iteration
        all_coords = numpy.concatenate(all_coords)
        with open(self.filename, 'w') as outfile:
            for i, frame in enumerate(all_coords):
                outfile.write("2\n")
                outfile.write("{0}\n".format(i))
                outfile.write("SOD {0:9.5f} {1:9.5f} {2:9.5f}\n".format(
                  float(frame[0,0]), float(frame[0,1]), float(frame[0,2])))
                outfile.write("CLA {0:9.5f} {1:9.5f} {2:9.5f}\n".format(
                  float(frame[1,0]), float(frame[1,1]), float(frame[1,2])))
        
        
        
            
