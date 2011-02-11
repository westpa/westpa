from __future__ import division; __metaclass__ = type

from collections import namedtuple

import numpy

InitialState = namedtuple('InitialState', ['label', 'initial_prob', 'recycle_prob', 'pcoord', 'bin'])
TargetState  = namedtuple('TargetState', ['label', 'bin'])

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
        
        # The initial distribution, a list of InitialDistribution namedtuples
        self.initial_states = None
        
        # Target states, a list of TargetState namedtuples
        self.target_states = None
        
        # Number of dimentions in progress coordinate data
        self.pcoord_ndim = 1
        
        # Length of progress coordinate data for each segment
        self.pcoord_len = 1
        
        # Data type of progress coordinate
        self.pcoord_dtype = numpy.float64
        
    def add_initial_state(self, label, initial_prob, recycle_prob, pcoord, bin = None):
        if self.initial_states is None:
            self.initial_states = []
        if bin is None:
            bin = self.region_set.get_bin_containing(pcoord)
        self.initial_states.append(InitialState(label, initial_prob, recycle_prob, pcoord, bin))
        
    def add_target_state(self, label, bin):
        try:
            self.target_states.append(TargetState(label, bin))
        except AttributeError:
            self.target_states = [TargetState(label,bin)]
        
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
