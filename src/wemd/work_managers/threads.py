from __future__ import division, print_function; __metaclass__ = type

import sys, logging, threading, multiprocessing
import Queue
import argparse
import wemd
from wemd.work_managers import WEMDWorkManager, WMFuture

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
            self.future._set_exception(e)
        else:
            self.future._set_result(result)
            
ShutdownSentinel = object()

class ThreadsWorkManager(WEMDWorkManager):
    '''A work manager using threads.'''
    
    def __init__(self):
        super(ThreadsWorkManager,self).__init__()
        self.n_workers = None
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
        
    def parse_aux_args(self, aux_args, do_help = False):
        parser = argparse.ArgumentParser(usage='%(prog)s [NON_WORK_MANAGER_OPTIONS] [OPTIONS]',
                                         add_help=False)
        
        runtime_config = wemd.rc.config
        group = parser.add_argument_group('threads work manager options')
                
        group.add_argument('-n', type=int, dest='n_workers', default=multiprocessing.cpu_count(),
                            help='Number of worker threads to run. (Default: %(default)s)')
        if do_help:
            parser.print_help()
            sys.exit(0)
        args, extra_args = parser.parse_known_args(aux_args)
        
        self.n_workers = runtime_config['work_manager.n_workers'] = args.n_workers
        
        return extra_args
        
    def startup(self):
        self.workers = [threading.Thread(target=self.runtask, args=[self.task_queue], name='worker-{:d}'.format(i)) 
                        for i in xrange(0, self.n_workers)]
        for thread in self.workers:
            log.debug('starting thread {!r}'.format(thread))
            thread.start()
        self.mode = self.MODE_MASTER
        return self.MODE_MASTER
        
    def shutdown(self, exit_code = 0):
        # Put one sentinel on the queue per worker, then wait for threads to terminate
        for i in xrange(0, self.n_workers):
            self.task_queue.put(ShutdownSentinel)
        for thread in self.workers:
            thread.join()
    
    