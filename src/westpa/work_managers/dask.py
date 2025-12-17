from itertools import islice

from dask import distributed

import westpa.work_managers as work_managers
from .core import WorkManager, WMFuture


class _FutureWrapper(WMFuture):

    def __init__(self, future: distributed.Future):
        super().__init__()
        self.future = future

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


class DaskWorkManager(WorkManager):
    """Submits computations to a Dask cluster.

    Parameters
    ----------
    client : dask.distributed.Client, optional
        Connection to a Dask cluster. Defaults to ``Client()``.

    """

    def __init__(self, client=None):
        super().__init__()
        self.client = client or distributed.Client()

    def submit(self, fn, args=None, kwargs=None):
        args = args or ()
        kwargs = kwargs or {}
        future = self.client.submit(fn, *args, **kwargs)
        return _FutureWrapper(future)

    def as_completed(self, futures):
        fmap = {future.to_dask(): future for future in futures}
        for future in distributed.as_completed(fmap.keys()):
            yield fmap[future]

    def submit_as_completed(self, task_generator, queue_size=None):
        futures = [self.submit(fn, args, kwargs) for (fn, args, kwargs) in islice(task_generator, queue_size)]
        pending = set(futures)
        while pending:
            completed, pending = distributed.wait(pending, return_when='FIRST_COMPLETED')
            futures = [self.submit(fn, args, kwargs) for (fn, args, kwargs) in islice(task_generator, len(completed))]
            pending.update(futures)
            for future in completed:
                yield future

    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        wmenv = wmenv or work_managers.environment.default_env
        group = parser.add_argument_group('options for Dask work manager')
        group.add_argument(
            wmenv.arg_flag('dask_scheduler_address'),
            metavar='SCHEDULER_ADDRESS',
            help="Address of the task scheduler (e.g., '127.0.0.1:8786').",
        )

    @classmethod
    def from_environ(cls, wmenv=None):
        wmenv = wmenv or work_managers.environment.default_env
        address = wmenv.get_val('dask_scheduler_address')
        client = distributed.Client(address)
        return cls(client)
