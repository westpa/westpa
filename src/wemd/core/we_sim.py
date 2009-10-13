from __future__ import division
import cPickle as pickle
import numpy
import time, datetime        
from wemd.core.particles import Particle, ParticleCollection
from wemd.core.segments import Segment
from wemd.core.errors import PropagationIncompleteError

ZERO_INTERVAL = datetime.timedelta(0)

__metaclass__ = type

import logging
log = logging.getLogger('wemd.core.we_sim')

class WESimIter:
    """
    Describes per-iteration information (summary or otherwise) for
    a WE simulation.
    """
    
    def __init__(self, we_iter = None, n_particles = None, norm = None,
                 cputime = None, walltime = None, data = None):
        self.we_iter = we_iter
        self.n_particles = n_particles
        self.norm = norm
        self.cputime = cputime
        self.walltime = walltime
        self.data = data
        
class WESimDriver:
    def __init__(self, runtime_config):
        self.we_driver = None
        self.work_manager = None
        self.runtime_config = runtime_config
        
        self.dbengine = runtime_config['data.db.engine']
        self.DBSession = runtime_config['data.db.sessionmaker']
        self.dbsession = None
        self.segments = None           
            
    def init_runtime(self):
        from string import Template
        self.runtime_config.require('data.state')
        drtemplate = self.runtime_config.setdefault('data.segrefs.template', 
                                                   'traj_segs/${we_iter}/${seg_id}')
        try:
            ctemplate = Template(drtemplate)
        except Exception, e:
            raise ConfigError('invalid data ref template', e)
        else:
            self.runtime_config['data.segrefs.ctemplate'] = ctemplate
        
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
        
    def make_data_ref(self, segment):
        template = self.runtime_config['data.segrefs.ctemplate']
        return template.substitute(segment.__dict__)
        
    def segments_to_particles(self, segments):
        """Convert (DB-oriented) segments into (WE-oriented) particles.
        Lineage information is not conserved, as it flows (one way) out of
        the WE routine."""
        particles = []
        for segment in segments:
            p = Particle(particle_id = segment.seg_id,
                         weight = segment.weight,
                         pcoord = segment.pcoord)
            particles.append(p)
        return particles
        
    def run_we(self):
        assert self.dbsession
        assert self.segments
        
        for segment in segments:
            if segment.status != Segment.SEG_STATUS_COMPLETE:
                raise PropagationIncompleteError('segment %s is not propagated'
                                                 % segment.seg_id)

        log.info('running WE on %d particles' % len(self.segments))
        current_particles = self.segments_to_particles(self.segments)
        
        # Perform actual WE calculation
        new_particles = self.we_driver.run_we(current_particles)
        new_we_iter = self.we_driver.current_iteration 
        
        # Convert particles to new DB segments
        new_segments = []
        self.dbsession.begin()
        sq = self.dbsession.query(Segment)
        for particle in new_particles:
            s = Segment(weight = particle.weight)
            s.we_iter = new_we_iter
            s.status = Segment.SEG_STATUS_PREPARED
            s.pcoord = None
            if particle.p_parent:
                s.p_parent = sq.get([particle.p_parent.particle_id])
                log.debug('segment %r parent is %r' % (s, s.p_parent))
            if particle.parents:
                s.parents = set(sq.filter(Segment.seg_id.in_([pp.particle_id for pp in particle.parents])).all())
                log.debug('segment %r parents are %r' % (s, s.parents))
            new_segments.append(s)
            self.dbsession.add(s)
        # Flush new segments to obtain segment IDs
        self.dbsession.flush()
        for segment in new_segments:
            segment.data_ref = self.make_data_ref(segment)
        # Record completed information about new segments
        self.dbsession.flush()
                        
        we_iter = WESimIter()
        we_iter.we_iter = new_we_iter
        we_iter.n_particles = len(new_segments)
        we_iter.norm = numpy.sum((seg.weight for seg in new_segments))
        # The "data" field is immutable, meaning that it will not get stored
        # unless a completely new object is specified for it
        we_data = {}
        we_data['bins_population'] = self.we_driver.bins_population
        we_data['bins_nparticles'] = self.we_driver.bins_nparticles
        if we_iter.we_iter > 0:
            we_data['bins_flux'] = self.we_driver.bins_flux
        we_iter.data = we_data
        
        self.dbsession.add(we_iter)
        self.dbsession.flush()
        self.dbsession.commit()
        self.segments = new_segments
             
    def init_sim(self, sim_config):
        import wemd.we_drivers
        from wemd.core import Segment, Particle
        from wemd.data_manager.schema import metadata
        
        for item in ('wemd.initial_particles', 'wemd.initial_pcoord'):
            sim_config.require(item)
            
        log.info('creating database tables')
        metadata.create_all(bind=self.dbengine)
        
        # Create and configure the WE driver
        self.we_driver = we_driver = wemd.we_drivers.make_we_driver(sim_config)
        we_driver.initialize(sim_config)
        
        # Create the initial segments
        log.info('creating initial segments')
        n_init = sim_config.get_int('wemd.initial_particles')
        pcoord = numpy.array([float(x) for x in 
                              sim_config.get_list('wemd.initial_pcoord')])        
        segments = [Segment(seg_id=i,
                            we_iter = 0, 
                            status = wemd.Segment.SEG_STATUS_COMPLETE,
                            weight=1.0/n_init,
                            pcoord = pcoord)
                    for i in xrange(1,n_init+1)]
        self.dbsession = self.DBSession()
        self.dbsession.begin()
        for segment in segments: 
            self.dbsession.add(segment)
        
        # Record dummy stats for the starting iteration
        stats = WESimIter()
        stats.we_iter = 0
        stats.cputime = stats.walltime = ZERO_INTERVAL
        stats.n_particles = len(segments)
        stats.norm = numpy.sum([seg.weight for seg in segments])
        self.dbsession.add(stats)
        
        # Record results to the database
        self.dbsession.commit()
        
        # Run one iteration of WE to assign particles to bins
        self.segments = segments
        self.run_we()
        self.dbsession.close()
        
    def run_sim(self):
        from sqlalchemy.orm import eagerload
        max_iterations = self.runtime_config.get_int('wemd.max_iterations')
        while self.we_driver.current_iteration <= max_iterations:
            log.info('WE iteration %d (of %d requested)'
                     % (self.we_driver.current_iteration, max_iterations))
            self.dbsession = self.DBSession()
            segs_remaining = self.dbsession.query(Segment)\
                .filter(Segment.we_iter == self.we_driver.current_iteration)\
                .filter(Segment.status != Segment.SEG_STATUS_COMPLETE)\
                .options(eagerload(Segment.p_parent))\
                .all()
            log.info('%d segments remaining in WE iteration %d' 
                     % (len(segs_remaining), self.we_driver.current_iteration))
            self.dbsession.begin()
            try:
                self.work_manager.propagate_segments(segs_remaining)
                self.dbsession.flush()
                self.dbsession.commit()
            except:
                self.dbsession.rollback()
                raise
            
            # check to see that all segments are SEG_STATUS_COMPLETE
            # or abort
            
            # load segments into self.segments
            # run we
            
            
            
