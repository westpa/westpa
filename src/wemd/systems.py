from __future__ import division; __metaclass__ = type

import numpy

class WEMDSystem:
    INITDIST_NAME = 0
    INITDIST_PROB = 1
    INITDIST_PCOORD = 2
    INITDIST_BIN = 3
    
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
        # The progress coordinate region set
        self.region_set = None
        
        # The initial distribution, a list of (name, probability, initial pcoord) tuples
        self.initial_distribution = None
        
        # Number of dimentions in progress coordinate data
        self.pcoord_ndim = 1
        
        # Length of progress coordinate data for each segment
        self.pcoord_len = 1
        
        # Data type of progress coordinate
        self.pcoord_dtype = numpy.float64
        
    def preprocess_segment(self, segment):
        '''Perform pre-processing on a given segment.  This is run by the worker immediately before propagation.'''
        pass
        
    def postprocess_segment(self, segment):
        '''Perform post-processing on a given segment.  This is run by the worker immediately after propagation.'''
        pass
