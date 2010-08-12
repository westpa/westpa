import cPickle as pickle
__metaclass__ = type

import logging
import time
import sys
log = logging.getLogger(__name__)

from wemd.rc import RC_SIM_CONFIG_KEY

class WESimManagerBase:
    """
    Simulation driver portions common to both master and workers
    """
    def __init__(self):
        # Runtime configuration; can be changed between invocations (to change file locations, etc)
        self.runtime_config = None

        # Static simulation configuration, set up and locked in at "wemdctl init"
        self.sim_config = None
        
        # The driver that actually propagates segments
        self.backend_driver = None
        
        self.data_manager = None
        self.we_driver = None
        
    def runtime_init(self, runtime_config, load_sim_config=True):
        self.runtime_config = runtime_config
        if load_sim_config:
            self.load_sim_config()
                    
    def sim_init(self, sim_config, sim_config_src):
        """Create the necessary state for a new simulation"""
        raise NotImplementedError
    
    def load_sim_config(self):
        """Load the static simulation configuration from disk"""
        self.runtime_config.require(RC_SIM_CONFIG_KEY)
        log.info("loading static simulation configuration from '%s'" % self.runtime_config[RC_SIM_CONFIG_KEY])
        self.sim_config = pickle.load(open(self.runtime_config[RC_SIM_CONFIG_KEY], 'rb'))
        
    def save_sim_config(self):
        """Save the static simulation configuration information to disk"""
        self.runtime_config.require(RC_SIM_CONFIG_KEY)
        log.info("saving static simulation configuration to '%s'" % self.runtime_config[RC_SIM_CONFIG_KEY])
        pickle.dump(self.sim_config, open(self.runtime_config[RC_SIM_CONFIG_KEY], 'wb'), pickle.HIGHEST_PROTOCOL)

    def load_data_manager(self):
        """Load and configure the data manager"""
        from wemd.data_manager import make_data_manager
        self.data_manager = make_data_manager(self.runtime_config)
    
    def load_we_driver(self):
        """Load and configure the WE driver"""
        from wemd.we_drivers import get_we_driver
        Driver = get_we_driver(self.sim_config['wemd.we_driver'])
        self.we_driver = Driver()
        self.we_driver.sim_config = self.sim_config
        self.we_driver.runtime_init(self.runtime_config)
        
    def load_backend_driver(self):
        from wemd.backend_drivers import get_backend_driver
        Driver = get_backend_driver(self.sim_config['backend.driver'])
        self.backend_driver = Driver()
        self.backend_driver.sim_config = self.sim_config
        self.backend_driver.runtime_init(self.runtime_config)
        
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
    
    def __init__(self):
        super(WESimMaster,self).__init__()
        self.we_iter = None
            
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
            
        if self.we_driver is None:
            self.load_we_driver()

        max_wallclock = self.max_wallclock   
        if( max_wallclock is not None):     
            we_cur_wallclock = time.time() - self.start_wallclock
            loop_start_time = loop_end_time = None
                    
        while self.continue_simulation():
            if( max_wallclock is not None):
                if( loop_end_time is not None):
                    loop_duration = loop_end_time - loop_start_time
                    we_cur_wallclock += loop_duration
                    if( we_cur_wallclock + loop_duration * 2.0 > max_wallclock ):
                        log.info('Shutdown so walltime does not exceed max wallclock:%r'%(max_wallclock))                        
                        self.shutdown(0)
                        sys.exit(0)

                loop_start_time = time.time()                        
                            
            self.prepare_iteration()
            self.backend_driver.pre_iter(self.we_iter)
            self.propagate_particles()
            self.run_we()
            self.backend_driver.post_iter(self.we_iter)
            self.finalize_iteration()
            
            if( max_wallclock is not None):
                loop_end_time = time.time()
                
    
    def save_sim_state(self):
        raise NotImplementedError
    
    def load_sim_state(self):
        raise NotImplementedError
 
class WESimClient(WESimManagerBase):
     def __init__(self):
         pass
     
     def load_sim_state(self):
         pass
    
from wemd.util.config_dict import ConfigError
def get_sim_manager(driver_name):
    from default import DefaultWEMaster
    if driver_name in ('', 'serial', 'default'):
        log.info('using default simulation manager')
        driver = DefaultWEMaster
    elif driver_name == 'mpi':
        log.info('using MPI simulation manager')
        from wemd.util import mpi as wemd_mpi
        
        wemd_mpi.init_mpi()
        
        from mpi import MPIWEMaster, MPIWEWorker
        if not wemd_mpi.is_mpi_active():
            log.warning('MPI environment not available; using serial driver')
            driver = DefaultWEMaster
        
        if wemd_mpi.is_rank_0():
            driver = MPIWEMaster
        else:
            driver = MPIWEWorker
    else:
        raise ConfigError('invalid simulation manager driver %r specified'
                          % driver_name)
    return driver
