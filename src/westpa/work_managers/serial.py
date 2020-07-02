import logging
import sys

from .core import WorkManager, WMFuture


log = logging.getLogger(__name__)


class SerialWorkManager(WorkManager):
    @classmethod
    def from_environ(cls, wmenv=None):
        return cls()

    def __init__(self):
        log.debug('initializing serial work manager')
        super().__init__()
        self.n_workers = 1

    def submit(self, fn, args=None, kwargs=None):
        ft = WMFuture()
        try:
            result = fn(*(args if args is not None else ()), **(kwargs if kwargs is not None else {}))
        except Exception as e:
            ft._set_exception(e, sys.exc_info()[2])
        else:
            ft._set_result(result)
        return ft
