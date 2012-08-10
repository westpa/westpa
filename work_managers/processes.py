from __future__ import division, print_function; __metaclass__ = type

import sys, logging, multiprocessing, threading, traceback, signal, os, random
import work_managers
from . import WorkManager, WMFuture

log = logging.getLogger(__name__)

# Tasks are tuples ('task', task_id, fn, args, kwargs).
# Results are tuples (rtype, task_id, payload) where rtype is 'result' or 'exception' and payload is the return value
# or exception, respectively.

task_shutdown_sentinel   = ('shutdown', None, None, (), {})
result_shutdown_sentinel = ('shutdown', None, None)

class ProcessWorkManager(WorkManager):
    '''A work manager using the ``multiprocessing`` module.'''
    
    @classmethod
    def from_environ(cls):
        return cls(work_managers.environment.get_worker_count())        
    
    def __init__(self, n_workers = None, shutdown_timeout = 1):
        super(ProcessWorkManager,self).__init__()
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self.workers = None
        self.task_queue = multiprocessing.Queue()
        self.result_queue = multiprocessing.Queue()
        self.receive_thread = None
        self.pending = None
        
        self.shutdown_received = False
        self.shutdown_timeout = shutdown_timeout or 1
        
    def task_loop(self):
        # Close standard input, so we don't get SIGINT from ^C
        try:
            sys.stdin.close()
        except Exception as e:
            log.info("can't close stdin: {}".format(e))

        # (re)initialize random number generator in this process
        random.seed()
        
        while not self.shutdown_received:
            message, task_id, fn, args, kwargs = self.task_queue.get()[:5]
            
            if message == 'shutdown':
                break
            
            try:
                result = fn(*args, **kwargs)
            except BaseException as e:
                result_tuple = ('exception', task_id, (e, traceback.format_exc()))
            else:
                result_tuple = ('result', task_id, result)
            self.result_queue.put(result_tuple)

        log.debug('exiting task_loop')
        return
        
    def results_loop(self):
        while not self.shutdown_received:
            message, task_id, payload = self.result_queue.get()[:3]
            
            if message == 'shutdown':
                break
            elif message == 'exception':
                future = self.pending.pop(task_id)
                future._set_exception(*payload)
            elif message == 'result':
                future = self.pending.pop(task_id)
                future._set_result(payload)
            else:
                raise AssertionError('unknown message {!r}'.format((message, task_id, payload)))

        log.debug('exiting results_loop')

    def submit(self, fn, *args, **kwargs):
        ft = WMFuture()
        log.debug('dispatching {!r}(*{!r},**{!r})'.format(fn,args,kwargs))
        self.pending[ft.task_id] = ft
        self.task_queue.put(('task', ft.task_id, fn, args, kwargs))        
        return ft
                
    def startup(self):
        self.workers = [multiprocessing.Process(target=self.task_loop, 
                                                name='worker-{:d}'.format(i)) for i in xrange(self.n_workers)]
        for worker in self.workers:
            worker.start()
            
        self.install_sigint_handler()
        self.pending = dict()

        self.receive_thread = threading.Thread(target=self.results_loop, name='receiver')
        self.receive_thread.daemon = True
        self.receive_thread.start()        
    
    def _empty_queues(self):
        while not self.task_queue.empty():
            try:
                self.task_queue.get(block=False)
            except multiprocessing.queues.Empty:
                break
            
        while not self.result_queue.empty():
            try:
                self.result_queue.get(block=False)
            except multiprocessing.queues.Empty:
                break        
        
    def shutdown(self):
        self._empty_queues()

        # Send shutdown signal
        for _i in xrange(self.n_workers):
            self.task_queue.put(task_shutdown_sentinel, block=False)
                    
        for worker in self.workers:
            worker.join(self.shutdown_timeout)
            if worker.is_alive():            
                log.debug('sending SIGINT to worker process {:d}'.format(worker.pid))
                os.kill(worker.pid, signal.SIGINT)
                worker.join(self.shutdown_timeout)
                if worker.is_alive():
                    log.warning('sending SIGKILL to worker process {:d}'.format(worker.pid))
                    os.kill(worker.pid, signal.SIGKILL)
                    worker.join()
                    
                log.debug('worker process {:d} terminated with code {:d}'.format(worker.pid, worker.exitcode))
            else:
                log.debug('worker process {:d} terminated gracefully with code {:d}'.format(worker.pid, worker.exitcode))
        
        self._empty_queues()
        self.result_queue.put(result_shutdown_sentinel)
                
    
    