from __future__ import division; __metaclass__ = type

import sys

import logging
log = logging.getLogger(__name__)

import cPickle as pickle

import wemd
from wemd.rc import RC_SIM_CONFIG_KEY, RC_SIM_STATE_KEY
from wemd.util import extloader

class WESimManager:
    """A state machine and communications broker"""
    def __init__(self, runtime_config = None):
        self.runtime_config = runtime_config or {}
        self.sim_config = {}
        
        self.data_manager = None
        self.we_driver = None
        self.work_manager = None
        self.pcoord_driver = None
        
        self.status_stream = sys.stdout
        
        
    def _sim_init(self, sim_config, sim_config_src):
        """Hook for derived sim managers to perform their own simulation initialization"""
        raise NotImplementedError
    
    def sim_init(self, sim_config, sim_config_src):
        """Initialize a new simulation using the configuration information in sim_config_src. Each
        dependent component (WE driver, etc.) processes information from sim_config_src and inserts the
        results in sim_config."""
        
        self.sim_config = sim_config = {}
        sim_config['__sim_manager__'] = self.__class__.__name__
        sim_config['__data_manager__'] = self.data_manager.__class__.__name__
        #sim_config['__we_driver__'] = self.we_driver.__class__.__name__
        
        self._sim_init(sim_config, sim_config_src)
        self.data_manager.sim_init(sim_config, sim_config_src)
        self.we_driver.sim_init(sim_config, sim_config_src)
        # work manager has no simulation state; it is completely run-time-configured
        
    def load_sim_state(self):
        """Load simulation state"""
        self.runtime_config.require(RC_SIM_STATE_KEY)
        state_file = self.runtime_config[RC_SIM_STATE_KEY]
        log.info('loading simulation state from %r' % state_file)
        state_data = pickle.load(open(state_file, 'rb'))
        
        for (key, dest) in [('data_manager', self.data_manager),
                            ('we_driver', self.we_driver),
                            ('work_manager', self.work_manager)]:
            dest.restore_state(state_data[key])
                
    def save_sim_state(self):
        """Save simulation state"""
        self.runtime_config.require(RC_SIM_STATE_KEY)
        state_file = self.runtime_config[RC_SIM_STATE_KEY]
        state_data = {}
        for (key, dest) in [('data_manager', self.data_manager),
                            ('we_driver', self.we_driver),
                            ('work_manager', self.work_manager)]:
            state_data[key] = dest.dump_state()
        
        log.info('saving simulation state to %r' % state_file)
        state_data = pickle.dump(state_data, open(state_file, 'wb'), pickle.HIGHEST_PROTOCOL)
    
    def load_plugin_component(self, callable_qualname, description, path = None):
        log.debug('loading %s from %s' % (description, callable_qualname))
        callable = extloader.get_callable(callable_qualname, path)
        component = callable(self)
        return component
        
    def load_data_manager(self):
        log.info('loading HDF5 data manager')
        self.data_manager = wemd.data_manager.WEMDDataManager(self)
    
    def load_we_driver(self):
        log.info('loading fixed boundary/fixed particle count WE driver')
        #from wemd.we_drivers.fixed_bins import FixedBinWEDriver
        #self.we_driver = FixedBinWEDriver(self)
        self.we_driver.runtime_init()
    
    def load_work_manager(self):
        self.work_manager = wemd.work_managers.default.DefaultWorkManager(self)
        self.work_manager.runtime_init()
    
    def run(self):
        """Begin (or continue) running a simulation"""
        raise NotImplementedError
    
    def shutdown(self, exit_code = 0):
        """Shut down a running simulation"""
        pass
