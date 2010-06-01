from __future__ import division
__metaclass__ = type
import cPickle as pickle
import numpy
import string, time, datetime
from copy import copy        

from wemd.sim_managers import WESimMaster
from wemd.core.we_sim import WESimIter
from wemd.core.particles import Particle, ParticleCollection
from wemd.core.segments import Segment
from wemd.core.errors import PropagationIncompleteError
from wemd.rc import RC_SIM_STATE_KEY 

import logging
log = logging.getLogger(__name__)

class DefaultWEMaster(WESimMaster):
    def __init__(self):
        super(DefaultWEMaster,self).__init__()
        self.max_iterations = 1
        self.worker_blocksize = 1
        
    def runtime_init(self, runtime_config, load_sim_config = True):
        super(DefaultWEMaster, self).runtime_init(runtime_config, load_sim_config)
        self.load_data_manager()
        self.max_iterations = runtime_config.get_int('limits.max_iterations', 1)
        self.worker_blocksize = runtime_config.get_int('backend.blocksize', 1)
                                                  
    def save_sim_state(self):
        state_filename = self.runtime_config[RC_SIM_STATE_KEY]
        log.info('saving state to %s' % state_filename)
        state_dict = {'we_driver': self.we_driver}
        log.debug('state info: %r' % state_dict)
        pickle.dump(state_dict, open(state_filename, 'wb'), -1)
    
    def load_sim_state(self):
        self.runtime_config.require(RC_SIM_STATE_KEY)
        state_filename = self.runtime_config[RC_SIM_STATE_KEY]
        log.info('loading state from %s' % state_filename)
        state_dict = pickle.load(open(state_filename))
        log.debug('state info: %r' % state_dict)
        self.we_driver = state_dict['we_driver']

    def sim_init(self, sim_config, sim_config_src):
        self.sim_config = sim_config        
        sim_config_src.require_all(('wemd.initial_particles', 'wemd.initial_pcoord',
                                    'backend.driver'))
        
        sim_config['wemd.we_driver'] = sim_config['bins.type'] = sim_config_src['bins.type'].lower()
        sim_config['backend.driver'] = sim_config_src['backend.driver'].lower()
        sim_config['wemd.initial_particles'] = sim_config_src.get_int('wemd.initial_particles')
        sim_config['wemd.initial_pcoord'] = numpy.array(sim_config_src.get_list('wemd.initial_pcoord', type=float))
            
        # Create the database
        self.data_manager.prepare_database()
        
        # Create and configure the backend driver
        self.load_backend_driver()
        self.backend_driver.sim_init(sim_config, sim_config_src)

        # Create and configure the WE driver
        self.load_we_driver()
        self.we_driver.sim_init(sim_config, sim_config_src)
        
        # The backend and WE driver are successfully initialized; save any static config to disk
        self.save_sim_config()
        
        
        # Create the initial segments
        log.info('creating initial segments')
        n_init = sim_config['wemd.initial_particles']
        pcoord_vals = sim_config['wemd.initial_pcoord']
        pcoord = numpy.empty((1,len(pcoord_vals)), numpy.float64)
        pcoord[0] = pcoord_vals        
        segments = [Segment(n_iter = 0, 
                            status = Segment.SEG_STATUS_COMPLETE,
                            weight=1.0/n_init,
                            pcoord = pcoord)
                    for i in xrange(1,n_init+1)]
        
        # Record dummy stats for the starting iteration
        self.we_iter = WESimIter()
        self.we_iter.binarray = self.we_driver.make_bins()
        self.we_iter.n_iter = 0
        self.we_iter.n_particles = len(segments)
        self.we_iter.norm = numpy.sum([seg.weight for seg in segments])
        self.data_manager.create_we_sim_iter(self.we_iter)
        
        # Run one iteration of WE to assign particles to bins (for bin data
        # only)
        self.run_we(segments)

        self.we_iter.data = {}
        self.we_iter.data['bin_boundaries'] = self.we_driver.bins.boundaries
        self.we_iter.data['bins_shape'] = self.we_driver.bins.shape
        self.we_iter.data['bin_ideal_num'] = self.we_driver.bins.ideal_num

        anparticles = self.we_driver.bins.nparticles_array()
        pops = self.we_driver.bins.population_array()
        self.we_iter.data['bin_n_particles'] = anparticles
        self.we_iter.data['bin_populations'] = pops
        self.data_manager.update_we_sim_iter(self.we_iter) 
                    
                    
    def run_we(self, initial_segments = None):
        current_iteration = self.we_driver.current_iteration
        
        # Get number of incomplete segments
        ninc = self.data_manager.num_incomplete_segments(self.we_iter)
        if ninc:
            raise PropagationIncompleteError('%d segments have not been completed'
                                             % ninc)
            
        # Get all completed segments
        if initial_segments:
            log.debug("Initial Segments")
            segments = initial_segments
        else:
            log.debug("Not Initial Segments")
            segments = self.data_manager.get_segments(Segment.n_iter == self.we_iter.n_iter,
                                                      load_p_parent = True)
        
        # Calculate WE iteration end time and accumulated CPU and wallclock time
        total_cputime = 0.0
        total_walltime = 0.0
        for segment in segments:
            total_cputime += segment.cputime or 0.0
            total_walltime += segment.walltime or 0.0
        self.we_iter.cputime = total_cputime
        self.we_iter.walltime = total_walltime
        
        we_starttime = time.clock() 
        log.info('running WE on %d particles' % len(segments))
        # Convert DB-oriented segments to WE-oriented particles
        current_particles = []
        for segment in segments:
            p = Particle(particle_id = segment.seg_id,
                         weight = segment.weight,
                         pcoord = segment.pcoord[-1,:])
            current_particles.append(p)
        current_particles = ParticleCollection(current_particles)

        norm = current_particles.norm
        log.info('norm = %.15g; error in norm %.6g' % (norm, norm-1))
                
        # Perform actual WE calculation
        we_int_starttime = time.clock()
        new_particles = self.we_driver.run_we(current_particles)
        we_int_stoptime = time.clock()
        new_we_iter = self.we_driver.current_iteration
        
        # Mark old segments as merged/recycled/continued
        if not initial_segments:
            segments_by_id = dict((segment.seg_id, segment) for segment in segments)
            
            for particle in self.we_driver.particles_merged:
                if segments_by_id.has_key(particle.particle_id):    
                    segments_by_id[particle.particle_id].endpoint_type = Segment.SEG_ENDPOINT_TYPE_MERGED
            for particle in self.we_driver.particles_escaped:
                segments_by_id[particle.particle_id].endpoint_type = Segment.SEG_ENDPOINT_TYPE_RECYCLED
                
            self.data_manager.update_segments(self.we_iter, segments)
        
        # Create storage for next WE iteration data        
        we_iter = WESimIter()
        we_iter.n_iter = new_we_iter
        
        # Convert particles (phase space points) to new propagation segments
        new_segments = []
        for particle in new_particles:
            s = Segment(n_iter = new_we_iter,
                        status = Segment.SEG_STATUS_PREPARED,
                        weight = particle.weight,
                        endpoint_type = Segment.SEG_ENDPOINT_TYPE_CONTINUATION,
                        pcoord = None)
            if current_iteration > 0:
                if particle.p_parent:

                    if particle.p_parent in [p for p in self.we_driver.particles_escaped]:
                        pass
                    else:
                        s.p_parent = self.data_manager.get_segment(particle.p_parent.particle_id)
                        log.debug('segment %r primary parent is %r' 
                                  % (s.seg_id or '(New)', s.p_parent.seg_id))
                else:
                    log.debug('segment %r has no primary parent; will restart in initial bin' % s)
                if particle.parents:
                    s.parents = set([self.data_manager.get_segment(pp.particle_id)
                                     for pp in particle.parents])
                    log.debug('segment %r parents are %r' 
                              % (s.seg_id or '(New)',
                                 [s2.particle_id for s2 in particle.parents]))
            new_segments.append(s)

        we_iter.n_particles = len(new_segments)
        we_iter.norm = numpy.sum((seg.weight for seg in new_segments))
        we_iter.binarray = self.we_driver.bins
        
        # Save the total probability that flowed off the edge of the world
        recycled_population = 0
        for particle in self.we_driver.particles_escaped:
            recycled_population += particle.weight 
        we_iter.data['recycled_population'] = recycled_population
        log.info('%.6g probability recycled' % recycled_population)
        log.info('bin populations:')
        log.info(str(self.we_driver.bins.population_array()))  

        
        self.data_manager.create_we_sim_iter(we_iter)
        self.data_manager.create_segments(we_iter, new_segments)
        we_endtime = time.clock()
        log.info('core WE procedures took %.2g seconds' % (we_int_stoptime - we_int_starttime))
        log.info('WE (including data management) took %.2g seconds'
                 % (we_endtime - we_starttime))
                     
    def continue_simulation(self):
        return bool(self.we_driver.current_iteration <= self.max_iterations)
    
    def prepare_iteration(self):
        self.we_iter = self.data_manager.get_we_sim_iter(self.we_driver.current_iteration)
        self.we_iter.starttime = datetime.datetime.now()
        self.we_iter.data = copy(self.we_iter.data)
        self.we_iter.data['bin_boundaries'] = self.we_driver.bins.boundaries
        self.we_iter.data['bins_shape'] = self.we_driver.bins.shape
        self.we_iter.data['bin_ideal_num'] = self.we_driver.bins.ideal_num

        anparticles = self.we_iter.data['bins_nparticles'] = self.we_driver.bins_nparticles
        self.we_iter.data['bins_population'] = self.we_driver.bins_population
        self.we_iter.data['bins_popchange'] = self.we_driver.bins_popchange 

        n_pop_bins = anparticles[anparticles != 0].size
        n_bins = len(self.we_driver.bins)
        
        log.info('%d / %d bins are populated' %( n_pop_bins, n_bins))

        self.data_manager.update_we_sim_iter(self.we_iter)
        
    def propagate_particles(self):
        current_iteration = self.we_iter.n_iter
        log.info('WE iteration %d (of %d requested)'
                 % (current_iteration, self.max_iterations))
        n_inc = self.data_manager.num_incomplete_segments(self.we_iter)
        log.info('%d segments remaining in WE iteration %d'
                 % (n_inc, current_iteration))
        #segments = self.data_manager.get_prepared_segments(self.we_iter)
        segments = self.data_manager.get_segments((Segment.n_iter == self.we_iter.n_iter)
                                                  &(Segment.status == Segment.SEG_STATUS_PREPARED),
                                                  load_p_parent = True)
        
        #for segment in segments: 
        for segment in map(None, *(iter(segments),) * self.worker_blocksize):
            if type(segment) is tuple:
                self.backend_driver.propagate_segments(list(segment))
                self.data_manager.update_segments(self.we_iter, list(segment))
            else:    
                self.backend_driver.propagate_segments([segment])
                self.data_manager.update_segments(self.we_iter, [segment])
        
    def finalize_iteration(self):
        self.we_iter.endtime = datetime.datetime.now()
        self.data_manager.update_we_sim_iter(self.we_iter)
        self.save_sim_state()
        
