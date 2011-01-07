import cPickle as pickle
__metaclass__ = type

import logging
import time
import sys
log = logging.getLogger(__name__)

from wemd.rc import RC_SIM_CONFIG_KEY
from wemd.util.config_dict import ConfigError
from wemd.rc import EX_ERROR
from wemd.core.segments import Segment

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
        
        # The driver that handles communication between workers (cores, nodes, etc)
        self.worker = None
        
        #hostname for tcp worker
        self.hostname = None
        self.cport = None
        self.key = None
        
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
        self.sim_config, self.sim_config_src = pickle.load(open(self.runtime_config[RC_SIM_CONFIG_KEY], 'rb'))
        
    def save_sim_config(self):
        """Save the static simulation configuration information to disk"""
        self.runtime_config.require(RC_SIM_CONFIG_KEY)
        log.info("saving static simulation configuration to '%s'" % self.runtime_config[RC_SIM_CONFIG_KEY])
        pickle.dump((self.sim_config,self.sim_config_src), open(self.runtime_config[RC_SIM_CONFIG_KEY], 'wb'), pickle.HIGHEST_PROTOCOL)

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
        
    def run(self):
        """Enter a running state, such as driving the simulation or waiting
        for network data to arrive for processing."""
        raise NotImplementedError

    def set_server_hostname(self, hostname):
        self.hostname = hostname
        
    def set_secret_key(self, key):
        self.key = key
        
    def set_cport(self, cport):
        self.cport = cport
        
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
        
    def run_we(self):
        """Bin, split, and merge particles."""
        raise NotImplementedError
    
    def continue_simulation(self):
        """Determine if another iteration will occur"""
        raise NotImplementedError
    
    def shutdown(self, exit_code=0):
        if self.data_manager is not None:
            self.data_manager.close_hdf5()
            
        if self.worker is not None:
            self.worker.shutdown(exit_code)
            
    def load_work_manager(self, driver_name):
        from wemd.work_managers.default import DefaultWorkManager
        if driver_name in ('', 'serial', 'default'):
            log.info('using default work manager')
            driver = DefaultWorkManager
        elif driver_name == 'mpi':
            log.info('using MPI work manager')
            from wemd.util import mpi as wemd_mpi
            
            wemd_mpi.init_mpi()
            
            from wemd.work_managers.mpi import MPIWEMaster, MPIWEWorker
            if not wemd_mpi.is_mpi_active():
                log.warning('MPI environment not available; using serial driver')
                driver = DefaultWEMaster
            
            if wemd_mpi.is_rank_0():
                driver = MPIWEMaster
            else:
                driver = MPIWEWorker
        elif driver_name == 'tcp_server':
            log.info('using TCP simulation manager (server)')
            from wemd.work_managers.tcpip import TCPWEMaster
            driver = TCPWEMaster
        elif driver_name == 'tcp_client':
            log.info('using TCP simulation manager (client)')
            from wemd.work_managers.tcpip import TCPWEWorker
            driver = TCPWEWorker
        else:
            raise ConfigError('invalid simulation manager driver %r specified'
                              % driver_name)
        return driver
        
    def run(self):
        
        if self.we_driver is None:
            self.load_we_driver()

        if self.worker is None:
            driver_name = self.runtime_config.get('work_manager.driver', 'serial').lower()
            WM = self.load_work_manager(driver_name)
            self.worker = WM(self.runtime_config, True)
            if self.hostname is not None:
                self.worker.set_server_hostname(self.hostname)
            if self.worker.worker_is_master() == False:
                if self.cport is not None:
                    self.worker.set_cport(self.cport)
            if self.key is not None:
                self.worker.set_secret_key(self.key)
                
        #only open db for write if worker is master
        if self.worker.worker_is_master() and self.data_manager is None:
            self.load_data_manager()
        
        self.worker.runtime_init(self.runtime_config, True)
    
        if self.worker.worker_is_master() == False:
            self.worker.run()
            return 
        
        max_wallclock = self.max_wallclock   
        if( max_wallclock is not None):     
            we_cur_wallclock = time.time() - self.start_wallclock
            loop_start_time = loop_end_time = None
        
        segments = None
        prev_max_seg_id = None
        max_seg_id = None
        incomplete_segs = False          
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
            
            current_iteration = self.we_iter.n_iter
            log.info('WE iteration %d (of %d requested)'
                     % (current_iteration, self.max_iterations))
            
            n_inc = self.data_manager.num_incomplete_segments(self.we_iter)
            log.info('%d segments remaining in WE iteration %d'
                     % (n_inc, current_iteration))
     
            if n_inc == 0: #can happen if a run crashed between iter
                segments = self.data_manager.get_segments(self.we_iter.n_iter,
                                  load_p_parent = True)
                max_seg_id = max([s.seg_id for s in segments])  
                segments = self.run_we(segments = segments)   
                self.worker.post_iter(self.we_iter)                
                self.finalize_iteration()
                continue
            
            if segments is None:
                segments = self.data_manager.get_segments(self.we_iter.n_iter, 
                                                          status_criteria = Segment.SEG_STATUS_PREPARED,
                                                          load_p_parent = True)
                if not segments: #crashed during run
                    #get incomplete segments
                    log.warn('restarting %d incomplete segments' % n_inc)
                    segments = self.data_manager.get_segments(self.we_iter.n_iter, 
                                                          status_criteria = Segment.SEG_STATUS_COMPLETE,
                                                          status_negate = True,
                                                          load_p_parent = True)
                    
                    max_seg_id = self.data_manager.get_max_seg_id(self.we_iter.n_iter)
                    incomplete_segs = True
                
            if incomplete_segs == False:
                prev_max_seg_id = max_seg_id
                seg_ids = [s.seg_id for s in segments if s.seg_id is not None]
                if len(seg_ids) > 0:
                    max_seg_id = max(seg_ids)
                else:
                    max_seg_id = None
            
            if max_seg_id is None:
                if prev_max_seg_id is None:
                    prev_max_seg_id = self.data_manager.get_max_seg_id(self.we_iter.n_iter - 1)

                for i in xrange(0,len(segments)):
                    segments[i].seg_id = prev_max_seg_id + i + 1
                max_seg_id = segments[-1].seg_id
                                                         
            segments = self.worker.propagate_particles(self.we_iter, segments)
                                         
            self.data_manager.update_segments(self.we_iter, list(segments))
                
            #the new segments are returned
            #the old ones are updated/used for we
            if incomplete_segs:
                #after run crashed, need to load all segments, not just incomplete ones
                #run_we will fail if the segments are still incomplete
                segments = self.data_manager.get_segments(self.we_iter.n_iter,
                                  load_p_parent = True)
                incomplete_segs = False
                
            segments = self.run_we(segments = segments)
            
            self.worker.post_iter(self.we_iter)                
            self.finalize_iteration()
                
            self.worker.finalize_iteration()
            
            #remove references to parents of segments
            #prevent memory problem
            for segment in segments:
                segment.parents = set()
                if segment.p_parent:
                        segment.p_parent.p_parent = None
                        
            if max_wallclock is not None:
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

def get_sim_manager(name):
    from default import DefaultSimManager
    return DefaultSimManager

