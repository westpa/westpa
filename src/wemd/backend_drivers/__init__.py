import logging
log = logging.getLogger(__name__)
__metaclass__ = type


class BackendDriver:
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config
        
    def propagate_segments(self, segments):
        raise NotImplementedError


from wemd.core import ConfigError
from executable import ExecutableBackend

def make_backend_driver(runtime_config):
    driver_name = runtime_config['backend.driver']
    if driver_name == 'executable':
        driver = ExecutableBackend(runtime_config)
    else:
        raise ConfigError('invalid backend driver (%s) specified' % driver_name)
    log.info('using %s propagation backend' % driver_name)
    return driver