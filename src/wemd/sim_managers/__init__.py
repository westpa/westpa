__metaclass__ = type

class WESimManagerBase:
    """
    Simulation driver portions common to both master and workers
    """
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config
        self.backend_driver = None
        self.we_driver = None
        
    def load_backend_driver(self):
        from wemd.backend_drivers import make_backend_driver
        self.backend_driver = make_backend_driver(self.runtime_config)
            
    def run(self):
        """Enter a running state, such as driving the simulation or waiting
        for network data to arrive for processing."""
        raise NotImplementedError
        
class WESimMaster(WESimManagerBase):
    """
    The overall driver of a WE simulation, responsible for all state
    I/O and overall coordination of workers.
    """    
    
    def load_we_driver(self, sim_config):
        from wemd.we_drivers import make_we_driver
        self.we_driver = make_we_driver(sim_config)

    def initialize_simulation(self, sim_config):
        """Create the necessary state for a new simulation"""
        raise NotImplementedError
            
    def prepare_iteration(self):
        """Prepare for the next WE iteration"""
        raise NotImplementedError
            
    def run_we(self):
        """Bin, split, and merge particles."""
        raise NotImplementedError
    
    def finalize_iteration(self):
        """Postprocess the propagated particles"""
        raise NotImplementedError
    
    def continue_simulation(self):
        """Determine if another iteration will occur"""
        raise NotImplementedError
    
    def run(self):
        while self.continue_simulation():
            self.prepare_iteration()
            self.propagate_particles()
            self.run_we()
            self.finalize_iteration()
    
    def save_state(self):
        raise NotImplementedError
    
    def restore_state(self):
        raise NotImplementedError

def make_sim_manager(runtime_config):
    from default import DefaultWEMaster
    return DefaultWEMaster(runtime_config)
