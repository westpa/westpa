__metaclass__ = type

import logging
log = logging.getLogger(__name__)

class WESimManagerBase:
    """
    Simulation driver portions common to both master and workers
    """
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config
        self.backend_driver = None
        self.backend_driver_name = None
        self.we_driver = None
        self.we_iter = None
        
    def load_backend_driver(self, sim_config=None):
        from wemd.backend_drivers import make_backend_driver
        self.backend_driver = make_backend_driver(self.backend_driver_name,
                                                  self.runtime_config, 
                                                  sim_config)
            
    def run(self):
        """Enter a running state, such as driving the simulation or waiting
        for network data to arrive for processing."""
        raise NotImplementedError

    def shutdown(self, exit_code=0):
        pass
        
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
        if self.backend_driver is None:
            self.load_backend_driver()
            
        while self.continue_simulation():
            self.prepare_iteration()
            self.backend_driver.pre_iter(self.we_iter)
            self.propagate_particles()
            self.run_we()
            self.backend_driver.post_iter(self.we_iter)
            self.finalize_iteration()
    
    def save_state(self):
        raise NotImplementedError
    
    def restore_state(self):
        raise NotImplementedError
    

from wemd.util.config_dict import ConfigError
def make_sim_manager(runtime_config):
    from default import DefaultWEMaster
    driver_name = runtime_config.get('sim_manager.driver', 'serial')
    driver_name = driver_name.lower()
    if driver_name in ('', 'serial', 'default'):
        log.info('using default simulation manager')
        return DefaultWEMaster(runtime_config)
    elif driver_name == 'mpi':
        log.info('using MPI simulation manager')
        from wemd.util import mpi as wemd_mpi
        
        wemd_mpi.init_mpi()
        
        from mpi import MPIWEMaster, MPIWEWorker
        if not wemd_mpi.is_mpi_active():
            log.warning('MPI environment not available; using serial driver')
            return DefaultWEMaster(runtime_config)
        
        if wemd_mpi.is_rank_0():
            return MPIWEMaster(runtime_config)
        else:
            return MPIWEWorker(runtime_config)
    else:
        raise ConfigError('invalid simulation manager driver %r specified'
                          % driver_name)