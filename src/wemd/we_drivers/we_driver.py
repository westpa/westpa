from __future__ import division
__metaclass__ = type

import re
from itertools import izip
import numpy
import math, random
from copy import copy

from wemd.core import ConfigError, ParticleExitError
from wemd.core.particles import Particle, ParticleCollection
from wemd.core.segments import Segment
from wemd.core.binarrays import BinArray

class WEDriver:
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        self.current_iteration = None
        self.segments = None
        self.segment_lineages = None
        self.bins = None
        self.bins_population = None
        self.bins_nparticles = None
        self.bins_flux = None
                    
    def make_bins(self):
        raise NotImplementedError
                                        
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
    
    def initialize(self, n_init, init_pcoord):
        self.current_iteration = 0
        pcoord = init_pcoord 
        segments = []
        for i in xrange(0, n_init):
            segment = Segment(we_iter = self.current_iteration,
                              seg_id = self.next_seg_id(),
                              weight = 1 / n_init,                              
                              status = Segment.SEG_STATUS_COMPLETE,
                              final_pcoord = pcoord
                              )
            segments.append(segment)
        self.segments = segments
            
    def run_we(self, segments):
        assert(segments)
        self.segments = None
        self.segment_lineages = None

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
            segment.data_ref = self.sim_manager.make_data_ref(segment)
            
            new_segments.append(segment)
            segment_lineages.extend((segment, parent)
                                    for parent in (particle.parents or [particle]))
            
        self.segments = new_segments
        self.segment_lineages = segment_lineages
            

class WESimIter:
    """
    Describes per-iteration information (summary or otherwise) for
    a WE simulation.
    """
    
    def __init__(self):
        self.we_iter = None
        self.n_particles = None
        self.norm = None
        self.cputime = None
        self.walltime = None
        self.data = {}
        
class WESimDriver:
    def __init__(self):
        self.we_info = None
        self.we_driver = None
    
    def initialize(self):
        pass
    
    def save_we_info(self):
        """Save per-iteration information
        """
        pass
    
    def load_we_info(self):
        """Load per-iteration information
        """
        pass
    
    def save_segments(self):
        pass

    def load_segments(self):
        pass
    
    def save_state(self):
        pass
    
    def load_state(self):
        pass

    def propagate_segments(self):
        pass
    
    def simulation_complete(self):
        pass
    
    def run_we(self):
        pass
    
    def run_sim(self):
        while not self.simulation_complete():
            self.work_manager.propagate_segments(self.current_segments)
            self.data_manager.save_segments(self.current_segments)
            self.we_driver.run_we()
            self.data_manager.save_segments(self.we_driver.segments)
            self.data_manager.save_lineage(self.we_driver.segment_lineages)
            self.state_manager.save_state(self.we_driver)
    
        
