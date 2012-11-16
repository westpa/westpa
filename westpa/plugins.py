from __future__ import division, print_function; __metaclass__ = type

class WESTPlugin:
    def __init__(self, rc, config):
        pass

    def validate_config(self):
        return True
    
    def register_sim_manager_hooks(self, sim_manager):
        pass
    
    
    """
    def register_data_manager_hooks(self, data_manager):
        pass
    
    def register_propagator_hooks(self, propagator):
        pass
    
    def register_system_driver_hooks(self, system_driver):
        pass
    
    def register_we_driver_hooks(self, we_driver):
        pass
    """

class PluginAcceptor:
    pass
    