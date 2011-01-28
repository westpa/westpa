from __future__ import division; __metaclass__ = type

import sys, os

import logging
log = logging.getLogger(__name__)

import cPickle as pickle

import wemd
from wemd.rc import RC_SIM_STATE_KEY
from wemd.util import extloader

class WESimManager:
    """A state machine and communications broker"""
    def __init__(self, runtime_config = None):
        self.runtime_config = runtime_config or {}
        
        self.data_manager = None
        self.we_driver = None
        self.work_manager = None
        
        self.propagator = None
        
        self.system = None
                
        self.status_stream = sys.stdout
                                
    def load_sim_state(self):
        """Load simulation state"""
        self.runtime_config.require(RC_SIM_STATE_KEY)
        state_file = self.runtime_config[RC_SIM_STATE_KEY]
        log.info('loading simulation state from %r' % state_file)
        state_data = pickle.load(open(state_file, 'rb'))
        
        for (key, dest) in [('data_manager', self.data_manager),
                            ('we_driver', self.we_driver),
                            ('work_manager', self.work_manager)]:
            try:
                dest.restore_state(state_data[key])
            except AttributeError as e:
                if 'NoneType' in str(e):
                    log.warning('cannot load %s state (%s not loaded)' % (key,key))
                else:
                    raise
                
    def save_sim_state(self):
        """Save simulation state"""
        self.runtime_config.require(RC_SIM_STATE_KEY)
        state_file = self.runtime_config[RC_SIM_STATE_KEY]
        state_data = {}
        for (key, dest) in [('data_manager', self.data_manager),
                            ('we_driver', self.we_driver),
                            ('work_manager', self.work_manager)]:

            try:
                state_data[key] = dest.dump_state()
            except AttributeError as e:
                if 'NoneType' in str(e):
                    log.warning('cannot save %s state (%s not loaded)' % (key,key))
                else:
                    raise
        
        log.info('saving simulation state to %r' % state_file)
        state_data = pickle.dump(state_data, open(state_file, 'wb'), pickle.HIGHEST_PROTOCOL)
            
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
        
    def load_system_driver(self):
        sysdrivername = self.runtime_config.require('system.system_driver')
        log.info('loading system driver %r' % sysdrivername)
        
        sysmodpath = self.runtime_config.get('system.module_path')
        if sysmodpath:
            pathinfo = [os.path.abspath(os.path.realpath(os.path.expanduser(os.path.expandvars(pathcomp))))
                        for pathcomp in sysmodpath.split(os.pathsep)]
        else:
            pathinfo = None
        
        self.system = extloader.get_object(sysdrivername, pathinfo)(self)
        log.debug('system driver is %r' % self.system)
                    
    def run(self):
        """Begin (or continue) running a simulation"""
        return self.work_manager.run()
    
    def shutdown(self, exit_code = 0):
        """Shut down a running simulation"""
        self.work_manager.shutdown(exit_code)

