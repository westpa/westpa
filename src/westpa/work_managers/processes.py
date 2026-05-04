import os
import sys
import random
import logging
import threading
import traceback
import multiprocessing
from multiprocessing.queues import Empty

import westpa.work_managers as work_managers
from .core import WorkManager, WMFuture

log = logging.getLogger(__name__)

# Tasks are tuples ('task', task_id, fn, args, kwargs).
# Results are tuples (rtype, task_id, payload) where rtype is 'result' or 'exception' and payload is the return value
# or exception, respectively.

task_shutdown_sentinel = ('shutdown', None, None, (), {})
result_shutdown_sentinel = ('shutdown', None, None)


class ProcessWorkManager(WorkManager):
    '''A work manager using the ``multiprocessing`` module.

    Notes
    -----

    On MacOS, as of Python 3.8 the default start method for multiprocessing launching new processes was changed from fork to spawn.
    On Linux, as of Python 3.14, the default start method for multiprocessing launching new processes was changed from fork to spawn.
    In general, spawn is more robust and efficient, however it requires serializability of everything being passed to the child process.
    In contrast, fork is much less memory efficient, as it makes a full copy of everything in the parent process.
    However, it does not require picklability.

    So, on MacOS and Linux, the method for launching new processes is explicitly changed to fork from the (UNIX-specific) default of spawn.
    Unix should default to fork.

    See https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods and
    https://docs.python.org/3/library/multiprocessing.html#the-spawn-and-forkserver-start-methods for more details.
    '''

    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env
        return cls(wmenv.get_val('n_workers', multiprocessing.cpu_count(), int))

    def __init__(self, n_workers=None, shutdown_timeout=1):
        super().__init__()

        try:
            if sys.platform in ['darwin', 'linux']:  # UNIX Platforms
                multiprocessing.set_start_method('fork')
                log.debug('setting multiprocessing start method to fork')
        except RuntimeError:
            log.debug('failed to set start method to fork')

        self.n_workers = n_workers or multiprocessing.cpu_count()
        self.workers = None
        self.task_queue = multiprocessing.Queue()
        self.result_queue = multiprocessing.Queue()
        self.receive_thread = None
        self.pending = None

        self.shutdown_received = threading.Event()
        self.shutdown_timeout = shutdown_timeout or 1

    def task_loop(self):
        # Close standard input, so we don't get SIGINT from ^C
        try:
            sys.stdin.close()
        except Exception as e:
            log.info("can't close stdin: {}".format(e))

        # (re)initialize random number generator in this process
        random.seed()

        while not self.shutdown_received.is_set():
            if not self.task_queue.empty():
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
        while not self.shutdown_received.is_set():
            if not self.result_queue.empty():
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
        return

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
        '''Empty self.task_queue and self.result_queue until queue is empty'''
        try:
            while True:
                self.task_queue.get_nowait()
        except (Empty, ValueError):
            log.debug('Empty task_queue')

        try:
            while True:
                self.result_queue.get_nowait()
        except (Empty, ValueError):
            log.debug('Empty result_queue')

    def shutdown(self):
        while self.running:
            log.debug('shutting down {!r}'.format(self))

            # Empty queues
            self._empty_queues()
            for _i in range(self.n_workers):
                self.task_queue.put_nowait(task_shutdown_sentinel)
            self.result_queue.put_nowait(result_shutdown_sentinel)

            # Signal shutdown Event to stop queue loops
            self.shutdown_received.set()
            self._empty_queues()

            # Terminating all workers
            for worker in self.workers:
                worker.join(self.shutdown_timeout)
                if worker.is_alive():
                    log.debug('sending SIGINT to worker process {:d}'.format(worker.pid))
                    worker.terminate()
                    worker.join(self.shutdown_timeout)
                    if worker.is_alive():
                        log.warning('sending SIGKILL to worker process {:d}'.format(worker.pid))
                        worker.kill()
                        worker.join()

                    log.debug('worker process {:d} terminated with code {:d}'.format(worker.pid, worker.exitcode))
                else:
                    log.debug('worker process {:d} terminated gracefully with code {:d}'.format(worker.pid, worker.exitcode))

                try:
                    worker.close()  # Release all resources
                except ValueError:
                    try:
                        if worker.is_alive():
                            log.debug('worker process {:d} could not be closed'.format(worker.pid))
                    except ValueError:
                        pass  # Already closed.

            self.running = False
            log.debug('Done shutting down the processes work manager')
