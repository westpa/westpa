from dask import distributed

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


class DaskWorkManager(WorkManager):
    """Submits computations to a Dask cluster.

    Parameters
    ----------
    client : dask.distributed.Client, optional

    """

    def __init__(self, client=None):
        super().__init__()
        self.client = client or distributed.Client()

    def submit(self, fn, args=None, kwargs=None):
        args = args or ()
        kwargs = kwargs or {}
        future = self.client.submit(fn, *args, **kwargs)
        return _FutureWrapper(future)
