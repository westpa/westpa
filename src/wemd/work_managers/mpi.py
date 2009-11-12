from __future__ import division
__metaclass__ = type
import os, sys, threading, Queue
from math import ceil
import cPickle as pickle
import logging
log = logging.getLogger(__name__)
from wemd.work_managers import WorkManager
from wemd.core.errors import PropagationIncompleteError
from mpi4py import MPI

# Ensure that the SQLAlchemy schema exists on this machine, even if a DB
# connection doesn't
import wemd.data_manager

class MPIWorkManager(WorkManager):
    def __init__(self):
        WorkManager.__init__(self)
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.rank
        self.nprocs = self.comm.size
    
    def initialize(self, backend_driver, runtime_config, sim_driver = None):
        super(MPIWorkManager,self).initialize(backend_driver, 
                                              runtime_config, 
                                              sim_driver)
        
    def propagate_segments(self, we_iter):
        if self.rank == 0:
            segments = self.sim_driver.q_incomplete_segments(we_iter).all()
        else:
            segments = []
            
        segments = self.comm.scatter(segments, root=0)
        
        log.info('processing %d segments' % len(segments))
        log.debug('segments: %s' % ([seg.seg_id for seg in segments],))
        for segment in segments:
            log.info('dispatching segment %s (weight %g)'
                     % (segment.seg_id, segment.weight))
            self.backend_driver.propagate_segment(segment)
            
        segments = self.comm.gather(segments, root=0)
        
        if self.rank == 0:
            dbsession = self.sim_driver.dbsession
            for segment in segments:
                dbsession.begin()
                try:
                    dbsession.merge(segment)
                    dbsession.add(segment)
                    dbsession.flush()
                except:
                    self.log.debug('error merging propagation results',
                                   exc_info = True)
                    dbsession.rollback()
                    raise
                else:
                    dbsession.commit()
        else:
            assert(not segments)
