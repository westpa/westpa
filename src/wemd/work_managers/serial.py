from __future__ import division; __metaclass__ = type

import logging

log = logging.getLogger(__name__)

from wemd.work_managers import WEMDWorkManager

class SerialWorkManager(WEMDWorkManager):
    def __init__(self):
        log.debug('initializing serial work manager')
        super(SerialWorkManager,self).__init__()
        
    def propagate(self, segments):
        propagator = self.propagator
        for segment in segments:
            propagator.propagate([segment])
