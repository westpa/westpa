from __future__ import division
__metaclass__ = type

import os, sys
import cPickle as pickle
import numpy
import logging
log = logging.getLogger(__name__)

from default import DefaultWEMaster
import wemd
import wemd.data_manager
import wemd.util.mpi
from wemd.sim_managers import WESimManagerBase

from mpi4py import MPI

class MPISimManager(WESimManagerBase):
    ROOT_TASK = 0
    
    TAG_WORKER_MESSAGE = 1
    MSG_EXIT = 1
    MSG_AWAIT_SEGMENT_SCATTER = 10
    
    def __init__(self, runtime_config):
        super(MPISimManager,self).__init__(runtime_config)
        self.comm = MPI.COMM_WORLD
                                
    def scatter_propagate_gather(self, segments):
        if segments:
            log.debug('scattering %d segments' % len(segments))
        elif segments is None:
            log.debug('awaiting scattered data')    
    
        segment = self.comm.scatter(segments, self.ROOT_TASK)
    
        if not segment:
            log.debug('no data to process')
        elif all(s == None for s in segment):
            log.debug('no data to process')
        else:
            for s in segment:
                log.debug('received segment %r via scatter' % s)
                
            if type(segment) is tuple:
                self.backend_driver.propagate_segments(list(segment))
            else:
                self.backend_driver.propagate_segments([segment])
                
            for s in segment:
                log.debug('dispatching segment %r via gather' % s)
        return self.comm.gather(segment, self.ROOT_TASK)
    
    def _get_message_name(self, msgno):
        msg_idx = dict((getattr(self, name), name) for name in dir(self)
                       if name.startswith('MSG_'))
        return msg_idx.get(msgno, str(msgno))
        
        
class MPIWEMaster(DefaultWEMaster, MPISimManager):
    def __init__(self, runtime_config):
        super(MPIWEMaster, self).__init__(runtime_config)
        
    def run(self):
        super(MPIWEMaster,self).run()
        self.send_directive(self.MSG_EXIT, 0)
        
    def send_directive(self, msg, data = None, block = True):
        pdata = pickle.dumps((msg,data), -1)
        requests = []
        if log.getEffectiveLevel() <= logging.DEBUG:
            msgname = self._get_message_name(msg)
        else:
            msgname = str(msg)
            
        for rank in xrange(0, self.comm.size):
            if rank == self.comm.rank: continue            
            # Send the data, and record the request handle
            requests.append(self.comm.Isend((pdata, len(pdata), MPI.BYTE), 
                                            rank,
                                            self.TAG_WORKER_MESSAGE))
            log.debug('dispatched message %s to rank %d' % (msgname, rank))
        if block:
            log.debug('waiting on delivery of all messages')
            MPI.Request.Waitall(requests)
            
        # pdata cannot be destroyed until actual send completes, so make sure
        # to hang onto a reference to it until *AFTER* the Waitall()
        # No way to do this with Python reference counting without using a
        # thread to dispatch the non-blocking sends and then wait on them,
        # then check whether the thread has terminated.  So, just pass out
        # a reference to pdata for the caller to hold onto if so desired.
        return requests, pdata
        
    def propagate_particles(self):
        current_iteration = self.we_driver.current_iteration
        log.info('WE iteration %d (of %d requested)'
                 % (current_iteration, self.max_iterations))
        n_inc = self.data_manager.num_incomplete_segments(current_iteration)
        n_workers = self.comm.size
        log.info('%d segments remaining in WE iteration %d'
                 % (n_inc, current_iteration))
        log.debug('dispatching to %d processes' % n_workers)
        
        prep_segments = self.data_manager.get_prepared_segments(current_iteration)
        while len(prep_segments) > 0:
            requests,pdata = self.send_directive(self.MSG_AWAIT_SEGMENT_SCATTER, 
                                                  block = False)
            #segments = prep_segments[0:n_workers]
            segments = prep_segments[0:n_workers*self.worker_blocksize]            
            if len(segments) < n_workers*self.worker_blocksize:
                segments = [None] * (n_workers*self.worker_blocksize - len(segments)) + segments
            log.debug('waiting on completed send of MSG_AWAIT_SEGMENT_SCATTER')                
            MPI.Request.Waitall(requests)            
            del pdata
            
            log.debug('scattering %d segments to %d workers' 
                      % (len(segments), n_workers))
                      
            # recast segment list into a list of tuples of len(n_workers)
            segments = map(None, *(iter(segments),) * self.worker_blocksize)
            
            segments = self.scatter_propagate_gather(segments)
            
            # flatten list of tuples into a single list
            segments = [j for i in segments for j in i]
            self.data_manager.update_segments(current_iteration, 
                                              [segment for segment in segments
                                               if segment is not None],
                                              )
            #del prep_segments[0:n_workers]
            del prep_segments[0:n_workers*self.worker_blocksize]
            
    def shutdown(self, exit_code=0):
        if exit_code != 0:
            log.error('terminating remote MPI processes')
        else:
            log.info('terminating remote MPI processes')
        self.send_directive(self.MSG_EXIT, exit_code, block = False)
            
class MPIWEWorker(MPISimManager):
    def __init__(self, runtime_config):
        super(MPIWEWorker, self).__init__(runtime_config)
        
    def restore_state(self):
        # No need for clients to restore state
        pass
    
    def unknown_message_abort(self, status, data):
        from wemd.environment import EX_COMM_ERROR
        log.fatal('unknown message %r (tag %d from rank %d) received'
                  % (data, status.tag, status.source))
        self.comm.Abort(EX_COMM_ERROR)
    
    def run(self):
        if self.backend_driver is None:
            self.load_backend_driver()
        log.debug('entering receive loop')
        while True:
            status = MPI.Status()
            data = self.comm.recv(source = self.ROOT_TASK,
                                  tag = MPI.ANY_TAG,
                                  status = status)
            (msgno, msgdata) = data        
            if status.tag == self.TAG_WORKER_MESSAGE:
                if log.getEffectiveLevel() <= logging.DEBUG:
                    log.debug('received message %s' % self._get_message_name(msgno))
                
                if msgno == self.MSG_AWAIT_SEGMENT_SCATTER:
                    log.debug('awaiting scatter')
                    self.scatter_propagate_gather(None)
                elif msgno == self.MSG_EXIT:
                    log.debug('exiting on directive from master process')
                    sys.exit(msgdata)
                else:
                    raise ValueError('invalid worker message %d' % msgno)
