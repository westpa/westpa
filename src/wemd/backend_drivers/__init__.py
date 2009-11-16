import logging
log = logging.getLogger(__name__)
__metaclass__ = type


class BackendDriver:
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config
        
    def pre_sim(self):
        pass
    
    def post_sim(self):
        pass
        
    def pre_iter(self, we_iter):
        pass
    
    def post_iter(self, we_iter):
        pass
    
    def pre_segment(self, segment):
        pass
    
    def post_segment(self, segment):
        pass
        
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