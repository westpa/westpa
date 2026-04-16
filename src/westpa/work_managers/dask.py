import logging
import os
from itertools import islice

import dask.distributed as distributed

import westpa
import westpa.work_managers as work_managers
from .core import WorkManager

log = logging.getLogger(__name__)


class _DaskFutureWrapper:
    # WMFuture-like interface to a dask.distributed.Future object.

    def __init__(self, future):
        super().__init__()
        self.future = future
        self._result = None

    def __repr__(self):
        return type(self).__name__ + '(' + repr(self.future) + ')'

    def __hash__(self):
        return hash(self.future)

    def get_result(self, discard=True):
        result = self.future.result()
        if not discard:
            self._result = result
        return result

    def get_exception(self):
        return self.future.exception()

    def get_traceback(self):
        return self.future.traceback()

    def is_done(self):
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


class DaskWorkManager(WorkManager):
    """Submits tasks to a Dask cluster.

    Parameters
    ----------
    client : dask.distributed.Client, optional
        Connection to a Dask cluster. If not provided, a ``LocalCluster``
        will be created.
    n_workers : int, optional
        Number of workers to use when creating a ``LocalCluster``.
        Ignored when connecting to an existing scheduler.

    """

    def __init__(self, client=None, n_workers=None):
        super().__init__()

        self.client = client
        self.supplied_n_workers = n_workers

    def startup(self):
        """
        Automatically called when entering a context manager.
        Usually called by each CLI tool.
        """
        if not self.running:
            if self.client is not None:
                self._local_cluster = None
                self.client = self.client
            else:
                self._local_cluster = distributed.LocalCluster(n_workers=self.supplied_n_workers)
                self.client = distributed.Client(self._local_cluster)
                log.info(f'Started local Dask cluster with {self.n_workers} workers')

            self.client.register_plugin(_ConfigSetter())
            self.running = True

    @property
    def n_workers(self):
        return len(self.client.scheduler_info()['workers'])

    def shutdown(self):
        """
        Automatically called when exiting context manager.
        """
        if self.running:
            if self._local_cluster is not None:
                self._local_cluster.close()
                self.client.shutdown()
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

    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        wmenv = wmenv or work_managers.environment.default_env
        group = parser.add_argument_group('options for Dask work manager')
        group.add_argument(
            wmenv.arg_flag('dask_scheduler_address'),
            metavar='SCHEDULER_ADDRESS',
            help="Address of a scheduler (e.g., '127.0.0.1:8786').",
        )
        group.add_argument(
            wmenv.arg_flag('dask_scheduler_file'),
            metavar='SCHEDULER_FILE',
            help="Path to a JSON file containing scheduler information.",
        )

    @classmethod
    def from_environ(cls, wmenv=None):
        # When no scheduler address or file is provided, a ``LocalCluster`` is
        # created automatically. The number of workers can be controlled with the
        # ``--n-workers`` CLI flag or the ``WM_N_WORKERS`` environment variable.
        wmenv = wmenv or work_managers.environment.default_env
        address = wmenv.get_val('dask_scheduler_address')
        scheduler_file = wmenv.get_val('dask_scheduler_file')

        n_workers_val = wmenv.get_val('n_workers')
        n_workers = int(n_workers_val) if n_workers_val is not None else None

        if address or scheduler_file:
            client = distributed.Client(address=address, scheduler_file=scheduler_file)
            return cls(client=client)
        else:
            return cls(n_workers=n_workers)
