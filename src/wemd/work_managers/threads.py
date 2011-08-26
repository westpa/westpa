from __future__ import division; __metaclass__ = type

import multiprocessing, threading, numpy, math, logging

log = logging.getLogger(__name__)

import wemd
from wemd.work_managers import WEMDWorkManager

# This is mostly for demonstration; serious parallelism probably needs processes, so that the
# global interpreter lock doesn't get in the way.

class ThreadedWorkManager(WEMDWorkManager):
    def __init__(self, propagator=None):
        log.debug('initializing threaded work manager')
        super(ThreadedWorkManager,self).__init__(propagator)
        self.cpu_count = multiprocessing.cpu_count()
        log.debug('cpu count: %d' % self.cpu_count)
        self.n_threads = wemd.rc.config.get_int('work_manager.n_threads', self.cpu_count)
        self.n_iter = None
        
        log.info('using %d threads for parallel propagation' % self.n_threads)
        
    def propagate(self, segments):
        n_rounds = int(math.ceil(len(segments) / self.n_threads))
        # Create an array large enough to hold all of our segments
        # Number of rows is n_threads, columns is n_rounds, large enough to
        # contain all of our segments.
        # Once that's done, each row contains the set of segments for the
        # corresponding thread to run
        segarray = numpy.empty((self.n_threads, n_rounds), numpy.object_)
        segarray.flat[0:len(segments)] = segments
        threads = [WorkerThread(self.propagator, segarray[ithread,:]) for ithread in xrange(0, self.n_threads)]
        for thread in threads:
            # Spawn threads, begin propagation in each
            thread.start()
        for thread in threads:
            # Wait on all threads 
            thread.join()            
        
class WorkerThread(threading.Thread):
    def __init__(self, propagator, segments):
        threading.Thread.__init__(self)
        self.propagator = propagator
        self.segments = segments
        
    def run(self):
        propagator = self.propagator
        segments = [segment for segment in self.segments if segment is not None]
        log.debug('propagating %d segment(s)' % len(segments))
        for segment in segments:
            assert segment is not None
            propagator.propagate([segment])
        log.debug('propagation complete')
