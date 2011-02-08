from __future__ import division; __metaclass__ = type

import numpy

class WEMDSystem:
    INITDIST_NAME = 0
    INITDIST_PROB = 1
    INITDIST_PCOORD = 2
    INITDIST_BIN = 3
    
    TARGET_NAME = 0
    TARGET_BIN = 1
    
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
        # The progress coordinate region set
        self.region_set = None
        
        # The initial distribution, a list of (name, probability, initial pcoord) tuples
        self.initial_states = None
        
        # Target states, a list of (name, bin) tuples
        self.target_states = None
        
        # Number of dimentions in progress coordinate data
        self.pcoord_ndim = 1
        
        # Length of progress coordinate data for each segment
        self.pcoord_len = 1
        
        # Data type of progress coordinate
        self.pcoord_dtype = numpy.float64
        
    def preprocess_iteration(self, n_iter, segments):
        '''Perform pre-processing on all segments for a new iteration.  This is
        run by the sim manager immediately prior to propagation.  Segment-based
        preprocessing (with preprocess_segments()) is to be preferred to this.'''
        pass
    
    def postprocess_iteration(self, n_iter, segments):
        '''Perform post-processing on all segments for an iteration.  This is
        run by the sim manager after propagation and prior to weighted
        ensemble.  Segment-based postprocessing (with postprocess_segments())
        is to be preferred to this.'''
        pass

    def preprocess_segments(self, segments):
        '''Perform pre-processing on a given set of segments.  This is run by a
        worker immediately before propagation of those segments.  Pre-processing
        of an entire iteration's set of segments at once should occur in
        preprocess_iteration().'''
        pass
        
    def postprocess_segments(self, segments):
        '''Perform post-processing on a given set of segments.  This is run by 
        a worker immediately after propagation of those segments.  Post-processing
        of an entire iteration's set of segments at once should occur in
        postprocess_iteration().'''
        pass
    
    def new_pcoord_array(self):
        return numpy.zeros((self.pcoord_len, self.pcoord_ndim), self.pcoord_dtype)
