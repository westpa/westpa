from __future__ import division
__metaclass__ = type

import os, sys
import cPickle as pickle
import numpy
import logging
log = logging.getLogger(__name__)

from default import DefaultWEMaster
import wemd
from wemd.sim_managers import WESimManagerBase

from mpi4py import MPI

class MPISimManager(WESimManagerBase):
    ROOT_TASK = 0
    
    TAG_EXIT = 1
    TAG_AWAIT_SCATTER = 10
    
    def __init__(self, runtime_config):
        super(MPISimManager,self).__init__(runtime_config)
        log.info('rank %d is %s:%d'
                 % (wemd.util.mpi.getrank(), wemd.util.mpi.getnodename(),
                    os.getpid()))
        self.load_backend_driver()
        self.comm = MPI.COMM_WORLD
        
    def recv_object(self, source=0, tag=0, status = None):
        rsize = numpy.array([0], 'i')
        log.debug('awaiting object size')
        self.comm.Recv((rsize, MPI.INT), source=source, tag=tag)
        log.debug('received object size %d' % rsize[0])
        log.debug('awaiting object')
        buf = numpy.array((rsize[0],), numpy.byte)
        self.comm.Recv((buf, MPI.BYTE), source=source, tag=tag)
        obj = pickle.loads(buf.data)
        return obj
        
    def scatter_propagate_gather(self, segments):
        if segments:
            log.debug('scattering %d segments' % len(segments))
        elif segments is None:
            log.debug('awaiting scattered data')
            
        segment = self.comm.scatter(segments, self.ROOT_TASK)
        if not segment:
            log.debug('no data to process')
        else:
            log.debug('received segment %r via scatter' % segment.seg_id)
            self.backend_driver.propagate_segments([segment])
            log.debug('dispatching segment %r via gather' % segment.seg_id)
        return self.comm.gather(segment, self.ROOT_TASK)
        
class MPIWEMaster(DefaultWEMaster, MPISimManager):
    def __init__(self, runtime_config):
        super(MPIWEMaster, self).__init__(runtime_config)
        
    def run(self):
        super(MPIWEMaster,self).run()
        self.send_directive(self.TAG_EXIT, 0)
        
    def send_directive(self, tag, data = None, block = True):
        size = numpy.array([0], 'i')
        bpdata = buffer(pickle.dumps(data))
        size[0] = len(bpdata)
        
        requests = []
        for rank in xrange(0, self.comm.size):
            if rank == self.comm.rank: continue
            # Send the size
            log.debug('sending object size %d to rank %d'  
                      % (size[0], rank))
            self.comm.Isend((size, MPI.INT), rank, tag)
            
            # Send the data, and record the request handle
            requests.append(self.comm.Isend((bpdata, MPI.BYTE), rank, tag))
            log.debug('dispatched message %d to rank %d' % (tag, rank))
        if block:
            log.debug('waiting on delivery of all messages')
            MPI.Request.Waitall(requests)
        else:
            return requests
        
    def propagate_particles(self):
        current_iteration = self.current_iteration
        log.info('WE iteration %d (of %d requested)'
                 % (current_iteration, self.max_iterations))
        q_inc = self.q_incomplete_segments(current_iteration)
        n_inc = q_inc.count()
        log.info('%d segments remaining in WE iteration %d'
                 % (n_inc, current_iteration))
        log.debug('dispatching to %d processes' % self.comm.size)
        
        # DANGER OF HANG HERE (if some segment can never finish...)
        while q_inc.count() > 0:
            requests = self.send_directive(self.TAG_AWAIT_SCATTER, block = False)
            segments = q_inc[0:self.comm.size]            
            if len(segments) < self.comm.size:
                segments = [None] * (self.comm.size - len(segments)) + segments
            log.debug('waiting on completed send of TAG_AWAIT_SCATTER')
            MPI.Request.Waitall(requests)
            
            log.debug('scattering %d segments to %d workers' 
                      % (len(segments), self.comm.size))
            segments = self.scatter_propagate_gather(segments)
            
            self.dbsession.begin()
            for segment in segments:
                if segment:
                    self.dbsession.merge(segment)
            try:
                self.dbsession.flush()
            except Exception, e:
                self.dbsession.rollback()
                raise e
            else:
                self.dbsession.commit()

class MPIWEWorker(MPISimManager):
    def __init__(self, runtime_config):
        super(MPIWEWorker, self).__init__(runtime_config)
        
    def restore_state(self):
        # No need for clients to restore state
        pass
        
    def run(self):
        log.info('entering receive loop')
        while True:
            status = MPI.Status()
            data = self.recv_object(source = self.ROOT_TASK, 
                                    tag=MPI.ANY_TAG, 
                                    status = status)
            if status.tag == self.TAG_AWAIT_SCATTER:
                log.debug('received scatter wait message')
                self.scatter_propagate_gather(None)
            elif status.tag == self.TAG_EXIT:
                log.debug('received exit message')
            else:
                log.fatal('unknown message (%d) received' % status.tag)
                self.comm.Abort(253)
