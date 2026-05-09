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

logger = logging.getLogger(__name__)


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
    client : dask.distributed.Client or dict[str, Any], optional
        Connection to a Dask cluster, or a dictionary of keyword arguments
        to be passed to ``Client``. If not provided, or if a dictionary
        without ``'address'`` or ``'scheduler_file'`` values is passed, a new
        ``LocalCluster`` will be created and managed by the work manager.
    **kwargs
        Keyword arguments to use when creating a ``LocalCluster``.
        Ignored if `client` specifies the cluster to use.

    """

    def __init__(self, client=None, **kwargs):
        super().__init__()

        self.client = client or {}
        self.kwargs = kwargs

        self._own_client = not isinstance(self.client, distributed.Client)
        self._local_cluster = None

    @property
    def n_workers(self):
        return len(self.client.scheduler_info()['workers'])

    def startup(self):
        if not self.running:
            if self._own_client:
                address = self.client.pop('address', None)
                scheduler_file = self.client.pop('scheduler_file', None)

                if address or scheduler_file:
                    # cluster created and managed by the user
                    self.client = distributed.Client(address=address, scheduler_file=scheduler_file, **self.client)
                else:
                    # cluster created and managed by the work manager
                    self._local_cluster = distributed.LocalCluster(**self.kwargs)
                    self.client = distributed.Client(self._local_cluster, **self.client)
                    logger.info(f'Started local Dask cluster with {self.n_workers} workers')

            self.client.register_plugin(_ConfigSetter(), name='config_setter')
            self.running = True

    def shutdown(self):
        if self.running:
            self.client.unregister_worker_plugin(name='config_setter')

            if self._own_client:
                self.client.close(timeout=5)

            if self._local_cluster is not None:
                self._local_cluster.close(timeout=5)
                self._local_cluster = None

            super().shutdown()

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

    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        wmenv = wmenv or work_managers.environment.default_env
        group = parser.add_argument_group('options for Dask work manager')
        group.add_argument(
            wmenv.arg_flag('dask_scheduler_address'),
            metavar='SCHEDULER_ADDRESS',
            help="Address of a Dask scheduler (e.g., '127.0.0.1:8786').",
        )
        group.add_argument(
            wmenv.arg_flag('dask_scheduler_file'),
            metavar='SCHEDULER_FILE',
            help="Path to a JSON file containing Dask scheduler information.",
        )
        group.add_argument(
            wmenv.arg_flag('dask_threads_per_worker'),
            metavar='THREADS_PER_WORKER',
            type=int,
            help="Number of threads per Dask worker. Ignored if SCHEDULER_ADDRESS or SCHEDULER_FILE is provided.",
        )
        group.add_argument(
            wmenv.arg_flag('dask_memory_limit'),
            metavar='MEMORY_LIMIT',
            type=str,
            help="Memory limit per Dask worker (e.g., '1GiB'). Ignored if SCHEDULER_ADDRESS or SCHEDULER_FILE is provided.",
        )

    @classmethod
    def from_environ(cls, wmenv=None):
        # When no scheduler address or file is provided, a ``LocalCluster`` is
        # created automatically. The number of workers can be controlled with the
        # ``--n-workers`` CLI flag or the ``WM_N_WORKERS`` environment variable.
        wmenv = wmenv or work_managers.environment.default_env
        client = {
            'address': wmenv.get_val('dask_scheduler_address'),
            'scheduler_file': wmenv.get_val('dask_scheduler_file'),
        }
        kwargs = {
            'n_workers': wmenv.get_val('n_workers', type_=int),
            'threads_per_worker': wmenv.get_val('dask_threads_per_worker', 1, type_=int),
            'memory_limit': wmenv.get_val('dask_memory_limit', 'auto'),
        }
        return cls(client, **kwargs)
