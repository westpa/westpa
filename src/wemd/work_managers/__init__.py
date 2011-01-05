import cPickle as pickle
__metaclass__ = type

import logging
import time
import sys
log = logging.getLogger(__name__)

from wemd.rc import RC_SIM_CONFIG_KEY 

class WEWorkManager:
    is_master = True
    is_worker = False

    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
                
        # The driver that actually propagates segments
        self.backend_driver = None
        
        self.we_iter = None
                
    def runtime_init(self):
        self.runtime_config = self.sim_manager.runtime_config
        self.sim_config = self.sim_manager.sim_config

        if self.backend_driver is None:
            self.load_backend_driver()            
                    
    def sim_init(self, sim_config, sim_config_src):
        """Create the necessary state for a new simulation"""
        raise NotImplementedError
        
    def load_backend_driver(self):
        from wemd.backend_drivers import get_backend_driver, get_backend_driver_by_file
        
        if( self.sim_config['backend.driver'] == 'file' ):
            Driver = get_backend_driver_by_file(self.runtime_config['backend.file.file'], 
                                                self.runtime_config['backend.file.class'])
        else:
            Driver = get_backend_driver(self.sim_config['backend.driver'])
            
        self.backend_driver = Driver()
        self.backend_driver.sim_config = self.sim_config
        self.backend_driver.sim_init(self.sim_config, self.sim_config_src)
        self.backend_driver.runtime_init(self.runtime_config)

    def propagate_particles(self, we_iter, segments):
        raise NotImplementedError
    
    def finalize_iteration(self):
        pass
    
    def prepare_iteration(self):
        pass
        
    def shutdown(self, exit_code=0):
        pass

class WEWorker:
    def __init__(self, sim_manager, work_manager):
        self.sim_manager = sim_manager
        self.work_manager = work_manager
        
    def get_status(self):
        """Return information about the current status of this worker"""
        raise NotImplementedError
        
    def propagate_particles(self, segments):
        """Use the backend to propagate particles"""
        raise NotImplementedError
    
    def prepare_iteration(self):
        """Set up worker state necessary for all segments to be run in a new iteration"""
        pass
    
    def finalize_iteration(self):
        """Finalize worker state after all segments in an iteration have been run"""
        pass
    
    def shutdown(self, exit_code = 0):
        """Gracefully terminate this worker"""
        pass
        
