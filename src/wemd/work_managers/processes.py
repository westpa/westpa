from __future__ import division; __metaclass__ = type

import multiprocessing, threading, numpy, math, logging

log = logging.getLogger(__name__)

from wemd.work_managers import WEMDWorkManager

# This is mostly for demonstration; serious parallelism probably needs processes, so that the
# global interpreter lock doesn't get in the way.

class ProcessWorkManager(WEMDWorkManager):
    def __init__(self, sim_manager):
        log.debug('initializing threaded work manager')
        super(ProcessWorkManager,self).__init__(sim_manager)
        self.cpu_count = multiprocessing.cpu_count()
        log.debug('cpu count: %d' % self.cpu_count)
        self.n_threads = sim_manager.runtime_config.get_int('work_manager.n_threads', self.cpu_count)
        self.n_iter = None
    
        log.info('using %d processes for parallel propagation' % self.n_threads)
        
    def propagate(self, segments):
        i = 0
        tasks = [[] for i in xrange(0, self.n_threads)]        
        for segment in map(None, *(iter(segments),) * self.n_threads):

            if type(segment) is tuple:
                tasks[i].append(segment)                
            else:    
                tasks[i].append((segment,))
                        
            i += 1
            i %= self.n_threads            

        #flatten list of tuples into a list
        processes = []
        segret = []
        queues = []
        for i in xrange(0, len(tasks)):
            seg = [j for k in tasks[i] for j in k if j is not None]
            if type(seg) is not list:
                seg = [seg]
            segret.append(seg)
            queue = multiprocessing.Queue()
            queues.append(queue)
            processes.extend([WorkerProcess(self.sim_manager, seg, queue)])
            
        for process in processes:
            # Spawn threads, begin propagation in each
            process.start()
        for process in processes:
            # Wait on all threads 
            process.join()      
            
        segret = []
        for i in xrange(0,len(queues)):
            segret.extend(queues[i].get())

        for i in xrange(0,len(segments)):
            segments[i] = segret[i]            
        
class WorkerProcess(multiprocessing.Process):
    def __init__(self, sim_manager, segments, queue):
        multiprocessing.Process.__init__(self)
        self.sim_manager = sim_manager
        self.segments = segments
        self.queue = queue
        
    def run(self):
        propagator = self.sim_manager.propagator
        system_driver = self.sim_manager.system
        
        log.debug('propagating %d segments' % len(self.segments))
        for i in xrange(0,len(self.segments)):
            segment = self.segments[i]
   
            if segment is None: continue

            propagator.propagate([segment])
            
            self.segments[i] = segment
            
        self.queue.put(self.segments)
        log.debug('propagation complete')
