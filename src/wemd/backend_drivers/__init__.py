import logging
log = logging.getLogger(__name__)
__metaclass__ = type


class BackendDriver:
    def __init__(self):
        self.runtime_config = None
        self.sim_config = None
    
    def sim_init(self, sim_config):
        pass
    
    def runtime_init(self, runtime_config):
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


def get_backend_driver(driver_name):
    from wemd.util.config_dict import ConfigError
    
    assert driver_name in ('executable', 'test', 'pyopenmm', 'pyopenmmmultiseg')
    log.info('using %s propagation backend' % driver_name)

    if driver_name == 'executable':
        from executable import ExecutableBackend
        return ExecutableBackend
    elif driver_name == 'test':
        from test import TestBackend
        return TestBackend
    elif driver_name == 'pyopenmm':
        try:
            from pyopenmm import PyOpenMMBackend
        except ImportError, e:
            raise ConfigError('pyopenmm backend driver unavailable (%s)' % e)
        else:
            return PyOpenMMBackend
    elif driver_name == 'pyopenmmmultiseg':
        try:
            from pyopenmm import PyOpenMMBackendMultiSeg
        except ImportError, e:
            raise ConfigError('pyopenmm-multiseg backend driver unavailable (%s)' % e)
        else:
            return PyOpenMMBackendMultiSeg
