from __future__ import division
import cPickle as pickle
import numpy
import time, datetime
from wemd.util import DBDataItem
from wemd.util.numpy_hacks import NumpyCmpSafeDict
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
        self.data = NumpyCmpSafeDict(data or dict())

class WESimDriver:
    def __init__(self):
        self.we_driver = None
        
        self.runtime_config = None
        self.dbengine = None
        self.DBSession = None
        self.dbsession = None
        self.segments = None           
            
    def init_runtime(self, runtime_config):
        from string import Template
        for item in ('data.state', 'data.db.url'):
            runtime_config.require(item)
            
        drtemplate = runtime_config.setdefault('data.segrefs.template', 
                                               'traj_segs/${we_iter}/${seg_id}')
        try:
            ctemplate = Template(drtemplate)
        except Exception, e:
            raise ConfigError('invalid data ref template', e)
        else:
            runtime_config['data.segrefs.ctemplate'] = ctemplate
            
        self.runtime_config = runtime_config
                
    def connect_db(self):
        import wemd.data_manager
        import sqlalchemy
        
        db_url = self.runtime_config['data.db.url']
        log.debug('connecting to %r' % db_url)
        self.dbengine = sqlalchemy.create_engine(db_url)
        self.DBSession = sqlalchemy.orm.sessionmaker(bind=self.dbengine,
                                                     autocommit=True,
                                                     autoflush=False)
        
    def save_state(self):
        state_dict = {'we_driver': self.we_driver}
        pickle.dump(state_dict, open(self.runtime_config['data.state'], 'wb'),
                    -1)
    
    def restore_state(self):
        state_dict = pickle.load(open(self.runtime_config['data.state']))
        self.sim_driver = state_dict['we_driver']
    
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
        we_iter.data['bins_population'] = self.we_driver.bins_population
        we_iter.data['bins_nparticles'] = self.we_driver.bins_nparticles
        if we_iter.we_iter > 0:
            we_iter.data['bins_flux'] = self.we_driver.bins_flux
        
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
            
        log.debug('creating database tables')
        self.dbengine.echo = True
        metadata.create_all(bind=self.dbengine)
        
        # Create and configure the WE driver
        self.we_driver = we_driver = wemd.we_drivers.make_we_driver(sim_config)
        we_driver.initialize(sim_config)
        
        # Create the initial segments
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
        
        
            
from wemd.core.particles import Particle, ParticleCollection
from wemd.core.segments import Segment
        
