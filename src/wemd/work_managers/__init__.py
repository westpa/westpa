import logging
log = logging.getLogger(__name__)
__metaclass__ = type

from wemd.core import Segment

class WorkManager:
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
    def propagate_segments(self, segments):
        raise NotImplementedError

from wemd.core import ConfigError

def make_work_manager(runtime_config):
    driver_name = runtime_config['work_manager.driver']
    if driver_name == 'serial':
        from serial import SerialWorkManager
        driver = SerialWorkManager()
    elif driver_name == 'mpi':
        from mpi import MPIWorkManager
        driver = MPIWorkManager()
    else:
        raise ConfigError('invalid work manager (%s) specified' % driver_name)
    log.info('using %s work manager' % driver_name)
    return driver

