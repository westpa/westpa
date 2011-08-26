from __future__ import division; __metaclass__ = type

import logging

log = logging.getLogger(__name__)

from wemd.work_managers import WEMDWorkManager

# This is mostly for demonstration; serious parallelism probably needs processes, so that the
# global interpreter lock doesn't get in the way.

class SerialWorkManager(WEMDWorkManager):
    def __init__(self, propagator=None):
        log.debug('initializing serial work manager')
        super(SerialWorkManager,self).__init__(propagator)
        
    def propagate(self, segments):
        propagator = self.propagator
        for segment in segments:
            propagator.propagate([segment])
