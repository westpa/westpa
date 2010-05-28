import logging
log = logging.getLogger(__name__)
__metaclass__ = type


class BackendDriver:
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config
                
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


def make_backend_driver(runtime_config):
    from wemd.core import ConfigError
    driver_name = runtime_config['backend.driver']
    if driver_name == 'executable':
        from executable import ExecutableBackend
        driver = ExecutableBackend(runtime_config)
    elif driver_name == 'test':
        from test import TestBackend
        driver = TestBackend(runtime_config)
    elif driver_name == 'pyopenmm':
        try:
            from pyopenmm import PyOpenMMBackend
        except ImportError, e:
            raise ConfigError('pyopenmm backend driver unavailable (%s)' % e)
        else:
            driver = PyOpenMMBackend(runtime_config)
    elif driver_name == 'pyopenmmmultiseg':
        try:
            from pyopenmm import PyOpenMMBackendMultiSeg
        except ImportError, e:
            raise ConfigError('pyopenmm-multiseg backend driver unavailable (%s)' % e)
        else:
            driver = PyOpenMMBackendMultiSeg(runtime_config)
    else:
        raise ConfigError('invalid backend driver (%s) specified' % driver_name)
    log.info('using %s propagation backend' % driver_name)
    return driver
