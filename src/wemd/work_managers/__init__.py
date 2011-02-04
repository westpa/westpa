import cPickle as pickle
__metaclass__ = type

import logging
import time
import sys
log = logging.getLogger(__name__)

class WEMDWorkManager:
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
        self.block_size = 1 # Number of segments this work manager can propagate in parallel
        # This is for information purposes only; the WEMD core will not dispatch to work managers
        # in chunks; self.block_size reports back to the core what chunking this work manager will *provide*
        
    def prepare_workers(self):
        '''Prepare workers for propagation'''
        pass
                        
    def prepare_iteration(self, n_iter, segments):
        self.n_iter = n_iter
        self.sim_manager.propagator.prepare_iteration(n_iter, segments)
                        
    def propagate(self, segments):
        raise NotImplementedError
    
    def finalize_iteration(self, n_iter, segments):
        self.sim_manager.propagator.finalize_iteration(n_iter, segments)
        
    def shutdown(self, exit_code=0):
        '''Cleanly shut down any active workers.'''
        pass

import serial, threads, processes
