from __future__ import division, print_function; __metaclass__ = type

import sys, logging, threading, multiprocessing
import Queue
from . import WorkManager, WMFuture
import work_managers

log = logging.getLogger(__name__)

class Task:        
    def __init__(self, fn, args, kwargs, future):
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.future = future
        
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as e:
            self.future._set_exception(e, sys.exc_info()[2])
        else:
            self.future._set_result(result)
            
ShutdownSentinel = object()

class ThreadsWorkManager(WorkManager):
    '''A work manager using threads.'''
    @classmethod
    def from_environ(cls): 
        return cls(work_managers.environment.get_worker_count())
    
    def __init__(self, n_workers = None):
        super(ThreadsWorkManager,self).__init__()
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self.workers = []
        self.task_queue = Queue.Queue()
        
    def runtask(self, task_queue):
        while True:
            task = task_queue.get()
            if task is ShutdownSentinel:
                return
            else:
                task.run()

    def submit(self, fn, *args, **kwargs):
        ft = WMFuture()
        task = Task(fn, args, kwargs, ft)
        self.task_queue.put(task)
        return ft
                
    def startup(self):
        self.workers = [threading.Thread(target=self.runtask, args=[self.task_queue], name='worker-{:d}'.format(i)) 
                        for i in xrange(0, self.n_workers)]
        for thread in self.workers:
            log.debug('starting thread {!r}'.format(thread))
            thread.start()
        
    def shutdown(self):
        # Put one sentinel on the queue per worker, then wait for threads to terminate
        for i in xrange(0, self.n_workers):
            self.task_queue.put(ShutdownSentinel)
        for thread in self.workers:
            thread.join()
    
    