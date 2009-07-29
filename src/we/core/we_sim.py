from __future__ import division
__metaclass__ = type

import re
from itertools import izip
import numpy
import math, random
from copy import copy

from we.core import ConfigError
from particles import Particle, ParticleCollection
from segments import Segment
from binarrays import BinArray

class WEError(Exception):
    pass

class PropagationIncompleteError(WEError):
    pass

class ParticleExitError(WEError):
    """Particle has exited the bin space
    args[0] = particle
    args[1] = tuple index into the bin limit array; None appears where the 
              particle did not fall into specified limits
              
    This may be raised as an error condition, or as a (rather inefficient)
    signal to process the particle in some way (after which the swarm will
    have to be re-binned from scratch)."""
    
    def __init__(self, particle, index):
        self.particle = particle
        self.index = index

    def __str__(self):
        return 'particle exited bin space; pcoord=%s, indices=%s' \
               % (self.particle.pcoord, self.index)  

class WESim:
    def __init__(self, data_manager = None):
        self.current_iteration = None
        self._next_seg_id = 0
        self.data_manager = data_manager
        self.bins = None
        self.bins_population = None
        self.bins_nparticles = None
        self.bins_flux = None
                
    def configure(self, sim_config):
        self.data_manager.config['we.initial_particles'] \
            = sim_config.get_int('we.initial_particles')
        self.data_manager.config['we.initial_pcoord'] \
            = numpy.array([float(val) for val in 
                           sim_config.get_list('we.initial_pcoord')]) 
        self.configure_data_refs(sim_config)
        self.configure_bin_params(sim_config)
        
    def configure_bin_params(self, sim_config):
        raise NotImplementedError
    
    def make_bins(self):
        raise NotImplementedError
        
    def configure_data_refs(self, config):
        from string import Template
        try:
            self.data_manager.config['data_refs.template'] \
                = Template(config['data_refs.template'])
        except Exception, e:
            raise ConfigError('could not compile data ref template', e)
        
    def find_next_seg_id(self, segments):
        self._next_seg_id = max(seg.seg_id for seg in segments) + 1
        return self._next_seg_id
        
    def next_seg_id(self, pretend = False):
        seg_id = self._next_seg_id
        if not pretend:
            self._next_seg_id += 1
        return seg_id 
    
    def make_data_ref(self, segment):
        dr_template = self.data_manager.config['data_refs.template']
        return dr_template.substitute(segment.__dict__)    
                    
    def distribute_particles(self, particles, binarray):
        boundaries = binarray.boundaries
        ndim = binarray.ndim
        for particle in particles:
            assert(len(particle.pcoord) == ndim)
            index = [None] * ndim
            for idim in xrange(0, ndim):
                for ibound in xrange(0, len(boundaries[idim])-1):
                    if boundaries[idim][ibound] <= particle.pcoord[idim] < boundaries[idim][ibound+1]:
                        index[idim] = ibound
            
            index = tuple(index)
            if None in index:
                self.on_particle_escape(particle, index)
                
            try:
                binarray.bins[index].add(particle)
            except IndexError:
                self.on_particle_escape(particle, index)
            
    def split_particles(self, particles, binarray):                
        for bin in binarray:
            bnorm = bin.norm
            ideal_weight = bnorm / bin.ideal_num
            for particle in list(bin):
                if particle.weight > bin.split_threshold*ideal_weight:
                    P = particle.__class__
                    m = int(math.ceil(particle.weight / ideal_weight))
                    assert(m >= 2)
                    
                    new_weight = particle.weight / m
                    new_particles = [P(seg_id = particle.seg_id,
                                       weight = new_weight,
                                       pcoord = copy(particle.pcoord),
                                       p_parent = particle,
                                       parents = [particle])]
                    new_particles.extend(P(seg_id = self.next_seg_id(),
                                           weight = new_weight,
                                           pcoord = copy(particle.pcoord),
                                           p_parent = particle,
                                           parents = [particle])
                                         for i in xrange(1, m))

                    self.on_particle_split(particle, new_particles)
                    
                    bin.discard(particle)
                    particles.discard(particle)
                    
                    bin.update(new_particles)
                    particles.update(new_particles)
    
    def merge_particles(self, particles, binarray):
        weight_key = (lambda p: p.weight)
        
        for bin in binarray:
            ideal_weight = bin.norm / bin.ideal_num
            min_weight = ideal_weight * bin.merge_threshold_min
            max_weight = ideal_weight * bin.merge_threshold_max
            
            mergable_particles = [particle for particle in bin
                                  if particle.weight <= min_weight]
            while len(mergable_particles) > 1:
                glom_src = []
                glom_weight = 0.0
                for particle in mergable_particles:
                    if glom_weight + particle.weight > max_weight:
                        break
                    else:
                        glom_src.append(particle)
                        glom_weight += particle.weight
                    
                if glom_src:
                    glom = Particle()                
                    u = random.uniform(0, glom_weight)
                    u_lower = glom_weight
                    u_upper = glom_weight
                    
                    glom_src = sorted(glom_src, key = weight_key, 
                                      reverse = True)
                    
                    for gparticle in glom_src:
                        u_lower -= gparticle.weight
                        if u_lower <= u < u_upper:
                            glom.pcoord = gparticle.pcoord
                            glom.p_parent = gparticle
                            break
                        else:
                            u_upper -= gparticle.weight
                    else:
                        # Should never happen, so throw a programmer-centric
                        # error
                        raise AssertionError('something wrong with '
                                             +'probabilities during merge')
                    
                    glom.weight = glom_weight
                    glom.parents = glom_src
                    glom.seg_id = glom.p_parent.seg_id
                    
                    self.on_particle_merge(glom_src, glom)
                    
                    bin.add(glom)
                    particles.add(glom)
                    
                    bin.difference_update(glom_src)
                    particles.difference_update(glom_src)
                # end if glom_src
                
                mergable_particles = [particle for particle in bin
                                      if particle.weight <= min_weight]
            # end while len(mergable_particles) > 1

    def on_particle_assign(self, particle, bin):
        pass

    def on_particle_escape(self, particle, index):
        raise ParticleExitError(particle, index)
    
    def on_particle_split(self, particle, children):
        pass
    
    def on_particle_merge(self, particles, glom):
        pass
    
    def initialize(self):
        self.current_iteration = 0
        self._next_seg_id = 0
        n_init = self.data_manager.config['we.initial_particles']
        pcoord = self.data_manager.config['we.initial_pcoord'] 
        segments = []
        for i in xrange(0, n_init):
            segment = Segment(we_iter = self.current_iteration,
                              seg_id = self.next_seg_id(),
                              weight = 1 / n_init,                              
                              status = Segment.SEG_STATUS_COMPLETE,
                              final_pcoord = pcoord
                              )
            segments.append(segment)
        self.data_manager.create_segments(segments)
        # Run an iteration of WE to perform initial bin assignments
        self.run_we()
            
    def run_we(self, pretend = False):
        if not self.data_manager.is_propagation_complete():
            raise PropagationIncompleteError()        
        segments = self.data_manager.get_segments(self.current_iteration)
        assert(segments)

        # Convert current segments to particles
        particles = ParticleCollection(Particle(seg_id = segment.seg_id,
                                                weight = segment.weight,
                                                pcoord = segment.final_pcoord)
                                       for segment in segments)

        # Bin particles
        last_bins = self.bins
        bins = self.bins = self.make_bins()
        # bin info stored here if time-dependent binning
        self.distribute_particles(particles, bins)
        
        # endpoint bin assignments stored here
        
        # Split/merge particles
        if self.current_iteration > 0:
            self.split_particles(particles, bins)
            self.merge_particles(particles, bins)
        
        # Calculate per-bin population and flux
        last_population = self.bins_population
        last_nparticles = self.bins_nparticles
        self.bins_population = bins.population_array()
        self.bins_nparticles = bins.nparticles_array()
        
        if self.current_iteration > 0:
            self.bins_flux = last_population - self.bins_population
            if not pretend:
                self.data_manager.record_data_item(self.current_iteration,
                                                   'bins_flux',
                                                   self.bins_flux)        
        if not pretend:
            self.data_manager.record_data_item(self.current_iteration,
                                               'bins_population',
                                               self.bins_population)

        # Convert particles to new segments
        if self.current_iteration > 0:
            self.find_next_seg_id(particles)
        self.current_iteration += 1
        new_segments = []
        segment_lineages = []
        for particle in particles:
            assert(particle.seg_id is not None)
            if not particle.p_parent:
                # Continuation
                p_parent_id = particle.seg_id
            else:
                p_parent_id = particle.p_parent.seg_id
            
            segment = Segment(we_iter = self.current_iteration,
                              seg_id = particle.seg_id,
                              p_parent_id = p_parent_id,
                              status = Segment.SEG_STATUS_PREPARED,
                              weight = particle.weight)
            segment.data_ref = self.make_data_ref(segment)
            
            new_segments.append(segment)
            segment_lineages.extend((segment, parent)
                                    for parent in (particle.parents or [particle]))
            
        if not pretend:
            self.data_manager.create_segments(new_segments)
            self.data_manager.record_lineage(segment_lineages)

