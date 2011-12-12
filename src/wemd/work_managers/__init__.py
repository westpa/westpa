import cPickle as pickle
__metaclass__ = type

import logging
import time
import sys
log = logging.getLogger(__name__)

import wemd

class WEMDWorkManager:
    def __init__(self, propagator=None):
        self.propagator = propagator
        self.shutdown_called = False
        self.mode = 'master'
                
    def parse_aux_args(self, aux_args, do_help = False):
        '''Parse any unprocessed command-line arguments, returning any arguments not proccessed
        by this object. By default, this does nothing except return the input.'''
        return aux_args
    
    def startup(self):
        '''Perform any necessary startup work, such as spawning clients'''
        pass
    
    def prepare_run(self):
        '''Prepare workers to distribute/receive work'''
        pass
                                
    def prepare_iteration(self, n_iter, segments):
        '''Prepare this work manager to run a new iteration. If overridden, this
        must call propagator.prepare_iteration().'''
        self.n_iter = n_iter
        self.propagator.prepare_iteration(n_iter, segments)
                        
    def propagate(self, segments):
        '''Propagate the given segments.'''
        raise NotImplementedError
    
    def finalize_iteration(self, n_iter, segments):
        '''Clean up at the end of an iteration.  If overridden, this must call
        propagator.finalize_iteration().'''
        self.propagator.finalize_iteration(n_iter, segments)
        
    def finalize_run(self):
        '''Clean up at the normal end of a run'''
        pass
        
    def shutdown(self, exit_code=0):
        '''Cleanly shut down any active workers.'''
        pass
    
    def run_worker(self):
        raise NotImplementedError('this work manager does not support dedicated workers; run as master')
    
import serial
