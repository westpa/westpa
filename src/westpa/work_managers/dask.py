import logging
import os
import sys
from itertools import islice

import westpa
import westpa.work_managers as work_managers
from .core import WorkManager

import dask

# This needs to be done before we import distributed, but after we import dask
if sys.platform in ['darwin', 'linux']:
    dask.config.set({'distributed.work.multiprocessing-method': 'fork'})
import dask.distributed as distributed

log = logging.getLogger(__name__)


class _DaskFutureWrapper:
    """WMFuture-like interface to a ``dask.distributed.Future`` object."""

    def __init__(self, future):
        super().__init__()
        self.future = future
        self._result = None

    def __repr__(self):
        return type(self).__name__ + '(' + repr(self.future) + ')'

    def __hash__(self):
        return hash(self.future)

    def get_result(self, discard=True):
        """Get result from ``distributed.client.Future``. By default,
        reference to future object will be removed after completion.
        """
        result = self.future.result()

        if not discard:
            self._result = result
        else:
            self.future = None
        return result

    def get_exception(self):
        return self.future.exception()

    def get_traceback(self):
        return self.future.traceback()

    def is_done(self):
        """Considered done if self.future reports done or no associated future object."""
        if self.future is None:
            return True
        else:
            return self.future.done()

    def wait(self):
        distributed.wait([self.future])

    @property
    def result(self):
        return self.get_result(discard=False)

    @property
    def exception(self):
        return self.get_exception()

    @property
    def traceback(self):
        return self.get_traceback()

    @property
    def done(self):
        return self.is_done()

    def to_dask(self):
        """Return ``dask.distributed.Future`` object"""
        return self.future


class _ConfigSetter(distributed.WorkerPlugin):
    # Distributes the client's WEST_SIM_ROOT environment variable and global
    # westpa.rc.config instance to workers.

    def __init__(self):
        self.sim_root = os.environ.get('WEST_SIM_ROOT') or os.getcwd()
        self.config = westpa.rc.config

    def setup(self, worker):
        os.environ['WEST_SIM_ROOT'] = self.sim_root
        westpa.rc.config = self.config

    def teardown(self, worker):
        del westpa.rc.config
        del self.config


class DaskWorkManager(WorkManager):
    """Submits tasks to a Dask cluster.

    Parameters
    ----------
    client : dask.distributed.Client or dict, optional
        Connection to a Dask cluster or a dictionary to the keyword arguments
        to be passed to ``dask.distributed.Client``. If not provided, a ``LocalCluster``
        will be created.
    n_workers : int, optional
        Number of workers to use when creating a ``LocalCluster``.
        Ignored when connecting to an existing scheduler.
    kwargs : dict, optional
        A dictionary with keyword arguments to be passed to ``LocalCluster``.
        ``kwargs`` should be of the following format::
            {'cluster': {'scheduler_file': '/path/to/file', 'address': 'localhost:12345'}}
    """

    def __init__(self, client=None, n_workers=None, threads_per_worker=None, **kwargs):
        super().__init__()

        self.client = client
        self.supplied_n_workers = n_workers
        self.supplied_n_threads = threads_per_worker or 1
        self.startup_kwargs = kwargs['cluster'] if 'cluster' in kwargs else kwargs

        print(f'{self.client=}')
        print(f'{self.startup_kwargs=}')

    def startup(self):
        """Automatically called when entering a context manager.
        Usually called by each CLI tool."""
        if not self.running:
            if self.client:
                if isinstance(self.client, dict):
                    print('a')
                    self.client = distributed.Client(
                        n_workers=self.supplied_n_workers, threads_per_worker=self.supplied_n_threads, **self.client
                    )
                    self._local_cluster = self.client.cluster
                else:
                    print('b')
                    self._local_cluster = distributed.LocalCluster(n_workers=self.supplied_n_workers, **self.startup_kwargs)
                    self.client = distributed.Client(self._local_cluster)
            else:
                print('c')
                self._local_cluster = distributed.LocalCluster(n_workers=self.supplied_n_workers, **self.startup_kwargs)
                self.client = distributed.Client(self._local_cluster)
                log.info(f'Started local Dask cluster with {self.n_workers} workers')

            self.client.register_plugin(_ConfigSetter(), name='config_setter')
            self.running = True

        if self._local_cluster is not None:
            print(self._local_cluster.workers)
            for nanny in self._local_cluster.workers.values():
                print(f'nanny: {nanny.pid} {repr(nanny)}')

        print(self.client.processing())
        print(self.client.scheduler_info())
        # print(self.client.dump_cluster_state(format='yaml'))

    @property
    def n_workers(self):
        return len(self.client.scheduler_info()['workers'])

    def shutdown(self):
        """Automatically called when exiting context manager."""
        if self.running:
            self.client.unregister_worker_plugin(name='config_setter')
            self.client.retire_workers(close_workers=True)
            self.client.scheduler.close()
            self.client.shutdown()

            if self._local_cluster is not None:
                for nanny in self._local_cluster.workers.values():
                    nanny.close(timeout=5, nanny=True)
                self._local_cluster.close()
                self._local_cluster = None

            super().shutdown()
            self.running = False

    def submit(self, fn, args=None, kwargs=None):
        args = args or ()
        kwargs = kwargs or {}
        future = self.client.submit(fn, *args, **kwargs, pure=False)
        return _DaskFutureWrapper(future)

    def as_completed(self, futures):
        fmap = {future.to_dask(): future for future in futures}
        for future in distributed.as_completed(fmap):
            yield fmap[future]

    def submit_as_completed(self, task_generator, queue_size=None):
        futures = [self.submit(fn, args, kwargs) for (fn, args, kwargs) in islice(task_generator, queue_size)]
        pending = {future.to_dask() for future in futures}
        while pending:
            completed, pending = distributed.wait(pending, return_when='FIRST_COMPLETED')
            futures = [self.submit(fn, args, kwargs) for (fn, args, kwargs) in islice(task_generator, len(completed))]
            pending |= {future.to_dask() for future in futures}
            for future in completed:
                yield _DaskFutureWrapper(future)

    def wait_any(self, futures):
        fmap = {future.to_dask(): future for future in futures}
        completed, pending = distributed.wait(list(fmap), return_when='FIRST_COMPLETED')
        return fmap[completed.pop()]

    def gather(self, futures):
        """Return all results associated to the current client for the given futures."""
        if isinstance(futures, distributed.Future):
            futures = [futures]

        fmap = {future.to_dask(): future for future in futures}
        return self.client.gather(fmap.keys())

    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        wmenv = wmenv or work_managers.environment.default_env
        group = parser.add_argument_group('options for Dask work manager')
        group.add_argument(
            wmenv.arg_flag('dask_scheduler_address'),
            metavar='DASK_SCHEDULER_ADDRESS',
            help="Address of a dask scheduler (e.g., '127.0.0.1:8786').",
        )
        group.add_argument(
            wmenv.arg_flag('dask_scheduler_file'),
            metavar='DASK_SCHEDULER_FILE',
            help="Path to a JSON file containing dask scheduler information.",
        )
        group.add_argument(
            wmenv.arg_flag('dask_n_threads_per_worker'),
            metavar='THREADS_PER_DASK_WORKER',
            type=int,
            help="Number of threads per dask worker.",
        )
        group.add_argument(
            wmenv.arg_flag('dask_memory_limit'),
            metavar='MEMORY_LIMIT',
            type=int,
            help="Memory limit per dask worker.",
        )

    @classmethod
    def from_environ(cls, wmenv=None):
        # When no scheduler address or file is provided, a ``LocalCluster`` is
        # created automatically. The number of workers can be controlled with the
        # ``--n-workers`` CLI flag or the ``WM_N_WORKERS`` environment variable.
        wmenv = wmenv or work_managers.environment.default_env

        # These are arguments passed to dask
        kwargs = {'client': {}, 'cluster': {}}
        kwargs['client']['address'] = wmenv.get_val('dask_scheduler_address')
        kwargs['client']['scheduler_file'] = wmenv.get_val('dask_scheduler_file')

        n_workers_val = wmenv.get_val('n_workers')
        kwargs['n_workers'] = int(n_workers_val) if n_workers_val is not None else None
        n_threads_val = wmenv.get_val('dask_n_threads_per_worker')
        kwargs['threads_per_worker'] = int(n_threads_val) if n_threads_val is not None else None
        memory_limit_val = wmenv.get_val('dask_memory_limit')
        kwargs['client']['memory_limit'] = memory_limit_val if memory_limit_val is not None else 'auto'

        return cls(**kwargs)
