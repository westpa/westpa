import logging
log = logging.getLogger(__name__)
__metaclass__ = type

class WorkManager:
    def __init__(self):
        self.backend_driver = None
        self.runtime_config = None
    
    def initialize(self, backend_driver, runtime_config):
        self.backend_driver = backend_driver
        self.runtime_config = runtime_config
    
    def propagate_segments(self, segments):
        raise NotImplementedError

from wemd.core import ConfigError

def make_work_manager(runtime_config):
    driver_name = runtime_config['work_manager.driver']
    if driver_name == 'serial':
        from serial import SerialWorkManager
        driver = SerialWorkManager()
    else:
        raise ConfigError('invalid work manager (%s) specified' % driver_name)
    log.info('using %s work manager' % driver_name)
    return driver

