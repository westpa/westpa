from __future__ import division; __metaclass__ = type

import logging, sys

log = logging.getLogger(__name__)

from wemd.work_managers import WEMDWorkManager, WMFuture

class SerialWorkManager(WEMDWorkManager):
    def __init__(self):
        log.debug('initializing serial work manager')
        super(SerialWorkManager,self).__init__()
        
    def submit(self, fn, *args, **kwargs):
        ft = WMFuture()
        try:
            result = fn(*args, **kwargs)
        except Exception as e:
            ft._set_exception(e, sys.exc_info()[2])
        else:
            ft._set_result(result)
        return ft
    