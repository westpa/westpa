import logging
import multiprocessing
import os
import random
import signal
import sys
import threading
import traceback

import westpa.work_managers as work_managers
from .core import WorkManager, WMFuture

log = logging.getLogger(__name__)

# Tasks are tuples ('task', task_id, fn, args, kwargs).
# Results are tuples (rtype, task_id, payload) where rtype is 'result' or 'exception' and payload is the return value
# or exception, respectively.

task_shutdown_sentinel = ('shutdown', None, None, (), {})
result_shutdown_sentinel = ('shutdown', None, None)


class ProcessWorkManager(WorkManager):
    '''A work manager using the ``multiprocessing`` module.'''

    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env
        return cls(wmenv.get_val('n_workers', multiprocessing.cpu_count(), int))

    def __init__(self, n_workers=None, shutdown_timeout=1):
        super().__init__()
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

    def submit(self, fn, args=None, kwargs=None):
        ft = WMFuture()
        log.debug('dispatching {!r}'.format(fn))
        self.pending[ft.task_id] = ft
        self.task_queue.put(('task', ft.task_id, fn, args or (), kwargs or {}))
        return ft

    def startup(self):
        from . import environment

        if not self.running:
            log.debug('starting up work manager {!r}'.format(self))
            self.running = True
            self.workers = [
                multiprocessing.Process(target=self.task_loop, name='worker-{:d}-{:x}'.format(i, id(self)))
                for i in range(self.n_workers)
            ]

            pi_name = '{}_PROCESS_INDEX'.format(environment.WMEnvironment.env_prefix)
            for iworker, worker in enumerate(self.workers):
                os.environ[pi_name] = str(iworker)
                worker.start()
            try:
                del os.environ[pi_name]
            except KeyError:
                pass

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
        if self.running:
            log.debug('shutting down {!r}'.format(self))
            self._empty_queues()

            # Send shutdown signal
            for _i in range(self.n_workers):
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
            self.running = False
