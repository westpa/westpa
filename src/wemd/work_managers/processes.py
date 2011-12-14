from __future__ import division, print_function; __metaclass__ = type

import sys, logging, multiprocessing, threading
import argparse
import wemd
from wemd.work_managers import WEMDWorkManager, WMFuture

log = logging.getLogger(__name__)

# Tasks are tuples ('task', task_id, fn, args, kwargs).
# Results are tuples (rtype, task_id, payload) where rtype is 'result' or 'exception' and payload is the return value
# or exception, respectively.

task_shutdown_sentinel   = ('shutdown', None, None, (), {})
result_shutdown_sentinel = ('shutdown', None, None)

class ProcessWorkManager(WEMDWorkManager):
    '''A work manager using the ``multiprocessing`` module.'''
    
    def __init__(self):
        super(ProcessWorkManager,self).__init__()
        self.n_workers = None
        self.workers = None
        self.task_queue = multiprocessing.Queue()
        self.result_queue = multiprocessing.Queue()
        self.receive_thread = None
        self.pending = None
        
        self.shutdown_timeout = 5
        
    def task_loop(self): 
        while True:
            message, task_id, fn, args, kwargs = self.task_queue.get()[:5]
            
            if message == 'shutdown':
                log.debug('shutting down worker loop')
                return
            
            try:
                result = fn(*args, **kwargs)
            except Exception as e:
                result_tuple = ('exception', task_id, e)
            else:
                result_tuple = ('result', task_id, result)
            self.result_queue.put(result_tuple)
        
    def results_loop(self):
        while True:
            message, task_id, payload = self.result_queue.get()[:3]
            
            if message == 'shutdown':
                log.debug('shutting down results collector')
                return
            elif message == 'exception':
                future = self.pending.pop(task_id)
                future._set_exception(payload)
            elif message == 'result':
                future = self.pending.pop(task_id)
                future._set_result(payload)
            else:
                raise AssertionError('unknown message {!r}'.format((message, task_id, payload)))

    def submit(self, fn, *args, **kwargs):
        ft = WMFuture()
        log.debug('dispatching {!r}(*{!r},**{!r})'.format(fn,args,kwargs))
        self.pending[ft.task_id] = ft
        self.task_queue.put(('task', ft.task_id, fn, args, kwargs))        
        return ft
        
    def parse_aux_args(self, aux_args, do_help = False):
        parser = argparse.ArgumentParser(usage='%(prog)s [NON_WORK_MANAGER_OPTIONS] [OPTIONS]',
                                         add_help=False)
        
        runtime_config = wemd.rc.config
        group = parser.add_argument_group('processes work manager options')
                
        group.add_argument('-n', type=int, dest='n_workers', default=multiprocessing.cpu_count(),
                            help='Number of worker processes to run. (Default: %(default)s)')
        if do_help:
            parser.print_help()
            sys.exit(0)
        args, extra_args = parser.parse_known_args(aux_args)
        
        self.n_workers = runtime_config['work_manager.n_workers'] = args.n_workers
        
        return extra_args
        
    def startup(self):
        self.workers = [multiprocessing.Process(target=self.task_loop, 
                                                name='worker-{:d}'.format(i)) for i in xrange(self.n_workers)]
        for worker in self.workers:
            worker.start()
            
        self.mode = self.MODE_MASTER
        self.pending = dict()

        self.receive_thread = threading.Thread(target=self.results_loop, name='receiver')
        self.receive_thread.daemon = True
        self.receive_thread.start()        
        return self.MODE_MASTER
        
    def shutdown(self, exit_code = 0):
        for _i in xrange(self.n_workers):
            self.task_queue.put(task_shutdown_sentinel, block=False)
        
        for worker in self.workers:
            log.debug('terminating process {:d}'.format(worker.pid))
            worker.join(self.shutdown_timeout)
            if worker.is_alive():
                log.warning('worker process {:d} did not exit; terminating'.format(worker.pid))
                worker.terminate()
            else:
                log.debug('worker process {:d} terminated gracefully with code {:d}'.format(worker.pid, worker.exitcode))
        
        self.result_queue.put(result_shutdown_sentinel)
                
    
    