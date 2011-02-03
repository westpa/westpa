from __future__ import division; __metaclass__ = type

import multiprocessing, threading, numpy, math, logging

log = logging.getLogger(__name__)

from wemd.work_managers import WEMDWorkManager

# This is mostly for demonstration; serious parallelism probably needs processes, so that the
# global interpreter lock doesn't get in the way.

class ThreadedWorkManager(WEMDWorkManager):
    def __init__(self, sim_manager):
        log.debug('initializing threaded work manager')
        super(ThreadedWorkManager,self).__init__(sim_manager)
        self.cpu_count = multiprocessing.cpu_count()
        log.debug('cpu count: %d' % self.cpu_count)
        self.n_threads = sim_manager.runtime_config.get_int('work_manager.n_threads', self.cpu_count)
        self.n_iter = None
        
        self.block_size = self.n_threads
        log.info('using %d threads for parallel propagation' % self.n_threads)
        
    def propagate_particles(self, segments):
        n_rounds = int(math.ceil(len(segments) / self.n_threads))
        segarray = numpy.empty((self.n_threads, n_rounds), numpy.object_)
        segarray.flat[0:len(segments)] = segments
        threads = [WorkerThread(self.sim_manager, segarray[ithread,:]) for ithread in xrange(0, self.n_threads)]
        for thread in threads:
            # Spawn threads, begin propagation in each
            thread.start()
        for thread in threads:
            # Wait on all threads 
            thread.join()            
        
class WorkerThread(threading.Thread):
    def __init__(self, sim_manager, segments):
        threading.Thread.__init__(self)
        self.sim_manager = sim_manager
        self.segments = segments
        
    def run(self):
        propagator = self.sim_manager.propagator
        system_driver = self.sim_manager.system
        log.debug('propagating %d segments' % len(self.segments))
        for segment in self.segments:
            if segment is None: continue
            log.debug('preparing segment %r'  % segment)
            propagator.prepare_segment(segment)
            log.debug('preprocessing segment %r' % segment)
            system_driver.preprocess_segment(segment)
            
            log.debug('propagating segment %r' % segment)
            propagator.propagate_segments([segment])
            
            log.debug('postprocessing segment %r' % segment)
            system_driver.postprocess_segment(segment)
            log.debug('finalizing segment %r' % segment)
            propagator.finalize_segment(segment)
        log.debug('propagation complete')
