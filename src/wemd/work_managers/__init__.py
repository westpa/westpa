import cPickle as pickle
__metaclass__ = type

import logging
import time
import sys
log = logging.getLogger(__name__)

from wemd.rc import RC_SIM_CONFIG_KEY 

class WEWorkManagerBase:
    """
    Simulation driver portions common to both master and workers
    """
    def __init__(self, runtime_config, load_sim_config = True):
        # Runtime configuration; can be changed between invocations (to change file locations, etc)
        self.runtime_config = None

        # Static simulation configuration, set up and locked in at "wemdctl init"
        self.sim_config = None
        self.sim_config_src = None
        
        # The driver that actually propagates segments
        self.backend_driver = None
        
        self.we_iter = None
                
    def runtime_init(self, runtime_config, load_sim_config=True):
        self.runtime_config = runtime_config

        #need sim_config to load backend driver
        if self.sim_config is None:
            self.load_sim_config()

        if self.backend_driver is None:
            self.load_backend_driver()            
                    
    def sim_init(self, sim_config, sim_config_src):
        """Create the necessary state for a new simulation"""
        raise NotImplementedError
    
    def load_sim_config(self):
        """Load the static simulation configuration from disk"""
        self.runtime_config.require(RC_SIM_CONFIG_KEY)
        log.info("loading static simulation configuration from '%s'" % self.runtime_config[RC_SIM_CONFIG_KEY])
        self.sim_config, self.sim_config_src = pickle.load(open(self.runtime_config[RC_SIM_CONFIG_KEY], 'rb'))

    def worker_is_master(self):
        return True        
    
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
        pass
    
    def finalize_iteration(self):
        """Runs right before next iter"""
        pass
    
    def prepare_iteration(self):
        pass
    
    def post_iter(self, we_iter):
        """Post iteration processing"""
        pass
    
    def shutdown(self, exit_code=0):
        pass
    