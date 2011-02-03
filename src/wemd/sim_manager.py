from __future__ import division; __metaclass__ = type

import sys, os, time, operator, numpy
from datetime import timedelta

import logging
log = logging.getLogger(__name__)

import cPickle as pickle

import wemd
from wemd.rc import RC_SIM_STATE_KEY
from wemd.util import extloader
from wemd.types import Segment, Particle

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
        
        self.state = {}
                                
    def load_sim_state(self):
        """Load simulation state"""
        self.runtime_config.require(RC_SIM_STATE_KEY)
        state_file = self.runtime_config[RC_SIM_STATE_KEY]
        log.info('loading simulation state from %r' % state_file)
        state_data = pickle.load(open(state_file, 'rb'))
        
        for (key, dest) in [('data_manager', self.data_manager),
                            ('we_driver', self.we_driver),]:
            try:
                dest.restore_state(state_data[key])
            except AttributeError as e:
                if 'NoneType' in str(e):
                    log.warning('cannot load %s state (%s not loaded)' % (key,key))
                else:
                    raise
        self.state.update(state_data['sim_manager'])
                
    def save_sim_state(self):
        """Save simulation state"""
        self.runtime_config.require(RC_SIM_STATE_KEY)
        state_file = self.runtime_config[RC_SIM_STATE_KEY]
        state_data = {}
        for (key, dest) in [('data_manager', self.data_manager),
                            ('we_driver', self.we_driver),]:

            try:
                state_data[key] = dest.dump_state()
            except AttributeError as e:
                if 'NoneType' in str(e):
                    log.warning('cannot save %s state (%s not loaded)' % (key,key))
                else:
                    raise
        state_data['sim_manager'] = self.state
        
        log.info('saving simulation state to %r' % state_file)
        state_data = pickle.dump(state_data, open(state_file, 'wb'), pickle.HIGHEST_PROTOCOL)
            
    def load_data_manager(self):
        log.info('loading HDF5 data manager')
        self.data_manager = wemd.data_manager.WEMDDataManager(self)
    
    def load_we_driver(self):
        log.info('loading default WE driver')
        self.we_driver = wemd.we_driver.WEMDWEDriver(self)
    
    def load_work_manager(self):
        self.work_manager = wemd.work_managers.threads.ThreadedWorkManager(self)
        
    def load_propagator(self):
        drivername = self.runtime_config.require('drivers.propagator')
        log.info('loading propagator %r' % drivername)
        pathinfo = self.runtime_config.get_pathlist('drivers.module_path')
        self.propagator = extloader.get_object(drivername, pathinfo)(self)
        log.debug('propagator is %r' % self.propagator)
        
    def load_system_driver(self):
        sysdrivername = self.runtime_config.require('system.system_driver')
        log.info('loading system driver %r' % sysdrivername)
        pathinfo = self.runtime_config.get_pathlist('system.module_path', None)        
        self.system = extloader.get_object(sysdrivername, pathinfo)(self)
        log.debug('system driver is %r' % self.system)
                    
    def run(self):
        """Begin (or continue) running a simulation"""
                
        # Have the work manager initialize
        #self.work_manager.prepare()
        
        # Set up internal timing
        run_start_time = time.time()
        max_walltime = self.runtime_config.get_interval('limits.max_wallclock', default=None, type=float)
        self.status_stream.write('Maximum wallclock time: %s\n' % timedelta(seconds=max_walltime or 0))
        last_iteration_runtime= 0
                        
        # Get segments
        n_iter = self.data_manager.current_iteration
        max_iter = self.runtime_config.get_int('limits.max_iterations', n_iter+1)
         
        
        
        # Guaranteed ordering by seg_id, so segments[seg_id] works for any valid seg_id for this iteration
        segments = sorted(self.data_manager.get_segments(n_iter = n_iter), key=operator.attrgetter('seg_id'))
        while n_iter <= max_iter: 
            self.status_stream.write('\n%s\n' % time.asctime())
            self.status_stream.write('Iteration %d (%d requested)\n' % (n_iter, max_iter))
            
            seg_weights = numpy.array(map(operator.attrgetter('weight'), segments))
            assert not (seg_weights == 0).any()
            norm = seg_weights.sum()
            self.status_stream.write('norm: %g, error in norm: %g\n' % (norm, norm-1))
            
            # Propagate segments
            segs_to_run = [segment for segment in segments if segment.status == Segment.SEG_STATUS_PREPARED]
            self.status_stream.write('%d of %d segments remain in iteration %d\n' % (len(segs_to_run), len(segments), n_iter))
            self.work_manager.prepare_iteration(n_iter)
            
            try:
                self.status_stream.flush()
            except AttributeError:
                pass
            
            self.work_manager.propagate_particles(segments)
            
            # Save results in case WE crashes
            self.data_manager.update_segments(n_iter, segments)
            
            # Check to ensure that all segments have been propagated
            failed_segments = [segment for segment in segments if segment.status != Segment.SEG_STATUS_COMPLETE]
            if failed_segments:
                self.status_stream.write('Propagation FAILED for %d segments:\n' % len(failed_segments))
                for failed_segment in failed_segments:
                    self.status_stream.write('  %d\n' % failed_segment.seg_id)
                raise RuntimeError('propagation failed for %d segments' % len(failed_segments))
            
            # Convert segments into particles representing their endpoints
            particles = [Particle(seg_id = segment.seg_id,
                                  weight = segment.weight,
                                  p_parent_id = None, # NOT the same as a segment's p_parent_id
                                  parent_ids = set(),
                                  pcoord = segment.pcoord[-1]) for segment in segments]
            
            # Run the weighted ensemble algorithm
            next_iter_particles = self.we_driver.run_we(particles, self.system.region_set)
                        
            # Mark old segments' endpoint types appropriately
            # Everything continues unless we are told otherwise
            n_recycled = 0
            P_recycled = 0.0
            n_merged = 0
            
            for segment in segments:
                segment.endpoint_type = Segment.SEG_ENDPOINT_TYPE_CONTINUES
                
            for seg_id in self.we_driver.recycle_terminations:
                segments[seg_id].endpoint_type = Segment.SEG_ENDPOINT_TYPE_RECYCLED
                n_recycled += 1
                P_recycled += segments[seg_id].weight
                
            for seg_id in self.we_driver.merge_terminations:
                segments[seg_id].endpoint_type = Segment.SEG_ENDPOINT_TYPE_MERGED
                n_merged += 1
                
            self.status_stream.write('%d particles (%g probability) recycled\n' % (n_recycled, P_recycled))
            self.status_stream.write('%d particles merged\n' % n_merged)
            self.status_stream.write('%d particles in next iteration\n' % len(next_iter_particles))
            
            self.data_manager.update_segments(n_iter, segments)
            self.data_manager.flush_backing()
            
            new_segments = []
            for particle in next_iter_particles:
                if log.isEnabledFor(logging.DEBUG):
                    log.debug('processing particle %r' % particle)
                segment = Segment(seg_id = None,
                                  weight = particle.weight,
                                  pcoord = numpy.expand_dims(particle.pcoord, 0),
                                  status = Segment.SEG_STATUS_PREPARED,
                                  )
                
                if particle.p_parent_id is None:
                    # Particle did not result from split or merge, but perhaps (if seg_id is negative) a recycle  
                    assert len(particle.parent_ids) == 0
                    assert particle.seg_id is not None
                    segment.p_parent_id = particle.seg_id
                    segment.parent_ids = set((particle.seg_id,))
                else:
                    # Particle did result from a split or a merge
                    assert len(particle.parent_ids) > 0
                    assert None not in particle.parent_ids
                    assert particle.seg_id is None
                    segment.p_parent_id = particle.p_parent_id
                    segment.parent_ids = set(particle.parent_ids)

                segment.n_parents = len(segment.parent_ids)
                new_segments.append(segment)
                    
            # Create new iteration group in HDF5            
            self.data_manager.prepare_iteration(n_iter+1, new_segments, 
                                                self.system.pcoord_ndim, self.system.pcoord_len, self.system.pcoord_dtype)
            
            n_iter += 1
            self.data_manager.current_iteration = n_iter
            segments = new_segments
            
            try:
                self.status_stream.flush()
            except AttributeError:
                pass
