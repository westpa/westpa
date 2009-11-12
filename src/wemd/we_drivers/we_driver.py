from __future__ import division
__metaclass__ = type

import logging
log = logging.getLogger('wemd.we_drivers.we_driver')

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
    def __init__(self, sim_config):
        self.current_iteration = 0
        self.particles = None
        self.particle_lineages = None
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
                    new_particles = [P(weight = new_weight,
                                       pcoord = copy(particle.pcoord),
                                       p_parent = particle,
                                       parents = [particle])]
                    new_particles.extend(P(weight = new_weight,
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
                    P = glom_src[0].__class__
                    glom = P()                
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
                    glom.particle_id = glom.p_parent.particle_id
                    
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
                
    def run_we(self, particles):
        assert(particles)

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
        
        # Assign parents as appropriate
        for particle in particles:
            if not particle.p_parent:
                # Continuation
                particle.p_parent = particle
                particle.parents = [particle]
        
        return particles

