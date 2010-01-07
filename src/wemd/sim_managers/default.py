from __future__ import division
__metaclass__ = type
import cPickle as pickle
import numpy
import string, time, datetime        

from wemd.sim_managers import WESimMaster
from wemd.core.we_sim import WESimIter
from wemd.core.particles import Particle, ParticleCollection
from wemd.core.segments import Segment
from wemd.core.errors import PropagationIncompleteError

import logging
log = logging.getLogger(__name__)

class DefaultWEMaster(WESimMaster):
    def __init__(self, runtime_config):
        super(DefaultWEMaster,self).__init__(runtime_config)
        if self.backend_driver is None: self.load_backend_driver()
        
        for key in (('data.state', 'backend.driver', 'data.storage_engine')):
            runtime_config.require(key)

        drtemplate = self.runtime_config.setdefault('data.segrefs.template', 
                                                   'traj_segs/${we_iter}/${seg_id}')

        ctemplate = string.Template(drtemplate)
        try:
            ctemplate.safe_substitute(dict())
        except ValueError, e:
            raise ConfigError('invalid data ref template %r' % drtemplate)
        else:
            self.runtime_config['data.segrefs.ctemplate'] = ctemplate

        self.max_iterations = runtime_config.get_int('limits.max_iterations', 1)
        # Eventually, add support for max_wallclock
        
        from wemd.data_managers import make_data_manager
        self.data_manager = make_data_manager(runtime_config)
        
    def make_data_ref(self, segment):
        template = self.runtime_config['data.segrefs.ctemplate']
        return template.safe_substitute(segment.__dict__)
                                  
    def save_state(self):
        state_filename = self.runtime_config['data.state']
        log.info('saving state to %s' % state_filename)
        state_dict = {'we_driver': self.we_driver}
        log.debug('state info: %r' % state_dict)
        pickle.dump(state_dict, open(state_filename, 'wb'), -1)
    
    def restore_state(self):
        state_filename = self.runtime_config['data.state']
        log.info('loading state from %s' % state_filename)
        state_dict = pickle.load(open(state_filename))
        log.debug('state info: %r' % state_dict)
        self.we_driver = state_dict['we_driver']

    def initialize_simulation(self, sim_config):
        import wemd.we_drivers
        from wemd.core import Segment, Particle
        
        for item in ('wemd.initial_particles', 'wemd.initial_pcoord'):
            sim_config.require(item)
            
        # Create the backing store
        self.data_manager.prepare_backing(sim_config)
                    
        # Create and configure the WE driver
        self.load_we_driver(sim_config)
        we_driver = self.we_driver
        
        # Create the initial segments
        log.info('creating initial segments')
        n_init = sim_config.get_int('wemd.initial_particles')
        pcoord_vals = [float(x) for x in 
                       sim_config.get_list('wemd.initial_pcoord')]
        pcoord = numpy.empty((1,len(pcoord_vals)), numpy.float64)
        pcoord[0] = pcoord_vals        
        segments = [Segment(n_iter = 0, 
                            status = wemd.Segment.SEG_STATUS_COMPLETE,
                            weight=1.0/n_init,
                            pcoord = pcoord)
                    for i in xrange(1,n_init+1)]
        
        # Record dummy stats for the starting iteration
        self.we_iter = WESimIter()
        self.we_iter.binarray = self.we_driver.make_bins()
        self.we_iter.n_iter = 0
        self.we_iter.n_particles = len(segments)
        self.we_iter.norm = numpy.sum([seg.weight for seg in segments])
        self.we_iter.segments = segments
        self.data_manager.create_we_sim_iter(self.we_iter)
        
        # Run one iteration of WE to assign particles to bins
        self.run_we()        
            
    def run_we(self):
        current_iteration = self.we_driver.current_iteration
        
        # Get number of incomplete segments
        ninc = self.data_manager.num_incomplete_segments(self.we_iter)
        if ninc:
            raise PropagationIncompleteError('%d segments have not been completed'
                                             % n_inc)
            
        # Get all completed segments
        segments = self.data_manager.get_segments(self.we_iter)
        
        # Calculate WE iteration end time and accumulated CPU and wallclock time
        total_cputime = 0.0
        total_walltime = 0.0
        for segment in segments:
            total_cputime += segment.cputime
            total_walltime += segment.walltime
        self.we_iter.cputime = total_cputime
        self.we_iter.walltime = total_walltime
                
        log.info('running WE on %d particles' % len(segments))
        # Convert DB-oriented segments to WE-oriented particles
        current_particles = []
        for segment in segments:
            p = Particle(particle_id = segment.seg_id,
                         weight = segment.weight,
                         pcoord = segment.pcoord)
            current_particles.append(p)
        current_particles = ParticleCollection(current_particles)

        norm = current_particles.norm
        log.info('norm = %.15g; error in norm %.6g' % (norm, norm-1))
                
        # Perform actual WE calculation
        new_particles = self.we_driver.run_we(current_particles)
        new_we_iter = self.we_driver.current_iteration 
        
        # Create storage for next WE iteration data        
        we_iter = WESimIter()
        we_iter.n_iter = new_we_iter
        
        # Convert particles to new propagation segments
        new_segments = []
        for particle in new_particles:
            s = Segment(weight = particle.weight)
            s.n_iter = new_we_iter
            s.status = Segment.SEG_STATUS_PREPARED
            s.pcoord = None
            if particle.p_parent:
                s.p_parent = Segment(n_iter = current_iteration,
                                     seg_id = particle.p_parent.particle_id)
                log.debug('segment %r primary parent is %r' 
                          % (s.seg_id or '(New)', s.p_parent.seg_id))
            if particle.parents:
                s.parents = set([Segment(n_iter = current_iteration,
                                         seg_id = pp.particle_id) for pp in particle.parents])
                log.debug('segment %r parents are %r' 
                          % (s.seg_id or '(New)',
                             [s2.particle_id for s2 in particle.parents]))
            new_segments.append(s)

        we_iter.segments = new_segments
        we_iter.n_particles = len(new_segments)
        we_iter.norm = numpy.sum((seg.weight for seg in new_segments))
        we_iter.binarray = self.we_driver.bins
        self.data_manager.create_we_sim_iter(we_iter)
                     
    def continue_simulation(self):
        return bool(self.we_driver.current_iteration <= self.max_iterations)
    
    def prepare_iteration(self):
        self.we_iter = self.data_manager.get_we_sim_iter(self.we_driver.current_iteration)
        
    def propagate_particles(self):
        current_iteration = self.we_iter.n_iter
        log.info('WE iteration %d (of %d requested)'
                 % (current_iteration, self.max_iterations))
        n_inc = self.data_manager.num_incomplete_segments(self.we_iter)
        log.info('%d segments remaining in WE iteration %d'
                 % (n_inc, current_iteration))
        for segment in self.data_manager.get_prepared_segments(self.we_iter):
            segments = [segment]
            self.backend_driver.propagate_segments(segments)
            
            self.dbsession.begin()
            try:
                self.dbsession.flush()
            except Exception, e:
                self.dbsession.rollback()
                raise
            else:
                self.dbsession.commit()
        
    def finalize_iteration(self):
        anparticles = self.we_driver.bins.nparticles_array()
        n_pop_bins = anparticles[anparticles != 0].size
        n_bins = len(self.we_driver.bins)
        log.info('%d / %d bins are populated' %( n_pop_bins, n_bins))
        self.save_state()
