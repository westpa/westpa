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
from sqlalchemy.orm import eagerload

ZERO_INTERVAL = datetime.timedelta(0)

import logging
log = logging.getLogger(__name__)

class DefaultWEMaster(WESimMaster):
    def __init__(self, runtime_config):
        super(DefaultWEMaster,self).__init__(runtime_config)
        if self.backend_driver is None: self.load_backend_driver()

        self.dbengine = None
        self.DBSession = None
        self._dbsession = None
        
        for key in (('data.state', 'data.db.url', 'backend.driver')):
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
            
        self.connect_db()
        
    def make_data_ref(self, segment):
        template = self.runtime_config['data.segrefs.ctemplate']
        return template.safe_substitute(segment.__dict__)
        
    current_iteration = property((lambda s: s.we_driver.current_iteration),
                                 None, None)
            
    def connect_db(self):
        db_url = self.runtime_config['data.db.url']
        log.info('connecting to %r' % db_url)
        
        import wemd.data_manager
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        
        self.dbengine = create_engine(db_url)
        self.DBSession = sessionmaker(bind=self.dbengine,
                                      autocommit = True,
                                      autoflush = False)
    
    def _get_dbsession(self):
        if self._dbsession is None:
            log.debug('implicitly instantiating new DB session')
            self._dbsession = self.DBSession()
        return self._dbsession
         
    def new_dbsession(self):
        if self._dbsession is not None:
            log.warn('implicitly closing open DB session to create new one')
            self.close_dbsession()
        return self._get_dbsession()
    
    def close_dbsession(self):
        self._dbsession.close()
        self._dbsession = None
        
    dbsession = property(_get_dbsession, None, close_dbsession, None)
         
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
            
    def q_incomplete_segments(self, we_iter = None):
        if we_iter is None: we_iter = self.current_iteration
        return self.dbsession.query(Segment)\
            .filter(Segment.we_iter == we_iter)\
            .filter(Segment.status != Segment.SEG_STATUS_COMPLETE)\
            .options(eagerload(Segment.p_parent))\
            .order_by(Segment.seg_id)
            
    def q_complete_segments(self, we_iter = None):
        if we_iter is None: we_iter = self.current_iteration
        return self.dbsession.query(Segment)\
                   .filter(Segment.we_iter == we_iter)\
                   .filter(Segment.status == Segment.SEG_STATUS_COMPLETE)
                                              
    def run_we(self):
        n_inc = self.q_incomplete_segments().count()
        if n_inc:
            raise PropagationIncompleteError('%d segments have not been completed'
                                             % n_inc)
        segments = self.q_complete_segments().all()
        
        # Record WE iteration end time and accumulated CPU and wallclock time
        import sqlalchemy
        SUM = sqlalchemy.func.sum
        
        self.dbsession.begin()
        sim_iter = self.dbsession.query(WESimIter).get([self.current_iteration])
        sim_iter.endtime = datetime.datetime.now()
        sim_iter.walltime = self.dbsession.query(SUM(Segment.walltime))\
                            .filter(Segment.we_iter == self.current_iteration)\
                            .scalar()
        sim_iter.cputime = self.dbsession.query(SUM(Segment.cputime))\
                           .filter(Segment.we_iter == self.current_iteration)\
                           .scalar()
        self.dbsession.flush()
        self.dbsession.commit()
        
        log.info('running WE on %d particles' % len(segments))
        # Convert DB-oriented segments to WE-oriented particles
        current_particles = []
        for segment in segments:
            p = Particle(particle_id = segment.seg_id,
                         weight = segment.weight,
                         pcoord = segment.pcoord)
            current_particles.append(p)
        current_particles = ParticleCollection(current_particles)

        log.info('norm = %.15g' % current_particles.norm)
        
        # Perform actual WE calculation
        new_particles = self.we_driver.run_we(current_particles)
        new_we_iter = self.we_driver.current_iteration 
        
        # Convert particles to new DB segments
        new_segments = []
        self.dbsession.begin()
        try:
            sq = self.dbsession.query(Segment)
            for particle in new_particles:
                s = Segment(weight = particle.weight)
                s.we_iter = new_we_iter
                s.status = Segment.SEG_STATUS_PREPARED
                s.pcoord = None
                if particle.p_parent:
                    s.p_parent = sq.get([particle.p_parent.particle_id])
                    log.debug('segment %r primary parent is %r' 
                              % (s.seg_id or '(New)', s.p_parent.seg_id))
                if particle.parents:
                    s.parents = set(sq.filter(Segment.seg_id.in_([pp.particle_id for pp in particle.parents])).all())
                    log.debug('segment %r parents are %r' 
                              % (s.seg_id or '(New)',
                                 [s2.particle_id for s2 in particle.parents]))
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
        except:
            log.debug('error in WE', exc_info = True)
            self.dbsession.rollback()
            raise
        else:
            self.dbsession.commit()
             
    def initialize_simulation(self, sim_config):
        import wemd.we_drivers
        from wemd.core import Segment, Particle
        from wemd.data_manager.schema import metadata
        
        for item in ('wemd.initial_particles', 'wemd.initial_pcoord'):
            sim_config.require(item)
            
        log.info('creating database tables')
        metadata.create_all(bind=self.dbengine)
        
        # Create and configure the WE driver
        self.load_we_driver(sim_config)
        we_driver = self.we_driver
        
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
        self.dbsession.begin()
        for segment in segments: 
            self.dbsession.add(segment)
        
        # Record dummy stats for the starting iteration
        stats = WESimIter()
        stats.we_iter = 0
        stats.cputime = stats.walltime = 0.0
        stats.n_particles = len(segments)
        stats.norm = numpy.sum([seg.weight for seg in segments])
        self.dbsession.add(stats)
        
        # Record results to the database
        self.dbsession.commit()
        
        # Run one iteration of WE to assign particles to bins
        self.run_we()        
        
    def continue_simulation(self):
        return bool(self.current_iteration <= self.max_iterations)
    
    def prepare_iteration(self):
        self.new_dbsession()
        sim_iter = self.dbsession.query(WESimIter).get([self.current_iteration])
        if sim_iter.starttime is None:
            sim_iter.starttime = datetime.datetime.now()
        self.dbsession.flush()
        
    def propagate_particles(self):
        current_iteration = self.current_iteration
        log.info('WE iteration %d (of %d requested)'
                 % (current_iteration, self.max_iterations))
        n_inc = self.q_incomplete_segments(current_iteration).count()
        log.info('%d segments remaining in WE iteration %d'
                 % (n_inc, current_iteration))
        for segment in self.q_incomplete_segments(current_iteration):
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
        self.save_state()
        self.close_dbsession()
