from __future__ import division, print_function; __metaclass__ = type

import sys, logging, multiprocessing, threading, traceback, signal, os
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
        
        self.shutdown_received = False
        self.shutdown_timeout = 5
            
    def child_handle_interrupt(self, signum, frame):
        # block sigint temporarily while we send it to child processes by sending it to
        # our own process group
        assert signum == signal.SIGINT
        self.shutdown_received = True
        prev_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        try:
            log.info('interrupted')
            pgid = os.getpgid(0)
            log.debug('sending SIGINT to process group {:d}'.format(pgid))
            # Kill child processes
            os.killpg(pgid, signal.SIGINT)
        finally:
            # Restore signal handler
            signal.signal(signal.SIGINT, prev_handler)
        
    def task_loop(self):
        try:
            sys.stdin.close()
        except Exception as e:
            log.info("can't close stdin: {}".format(e))
            
        signal.signal(signal.SIGINT, self.child_handle_interrupt)
        # Become our own process group
        os.setpgid(0,0) 
        while not self.shutdown_received:
            message, task_id, fn, args, kwargs = self.task_queue.get()[:5]
            
            if message == 'shutdown':
                break
            
            try:
                result = fn(*args, **kwargs)
            except Exception as e:
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
        
    def shutdown(self, exit_code = 0):
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
                
    
    