from __future__ import division
__metaclass__ = type

import logging
log = logging.getLogger('wemd.we_drivers.we_driver')

import math, random
from copy import copy

class WEDriver:
    def __init__(self, sim_config):
        self.current_iteration = 0
        
        # Newly-created particles (for whatever reason)
        self.particles_created = None
        
        # Disappearing particles
        self.particles_split = None
        self.particles_merged = None
        self.particles_escaped = None
        
        self.bins = None
        self.bins_population = None
        self.bins_nparticles = None
        self.bins_popchange = None
                            
    def make_bins(self):
        """Create an array of ParticleCollection objects appropriate to this WE
        method."""
        raise NotImplementedError
                                        
    def distribute_particles(self, particles, binarray):
        """Assign particles to bins.  Particles outside of the bin space must
        be assigned to the `particles_escaped` collection and may optionally
        have a `bin_index` attribute added to them, which will indicate which
        bin dimensions contain the particles fall into and which dimensions
        cannot contain the particles.
        """
         
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
                particle.bin_index = index
                self.particles_escaped.add(particle)
            else:                
                try:
                    binarray.bins[index].add(particle)
                except IndexError:
                    particle.bin_index = index
                    self.particles_escaped.add(particle)
            
    def split_particles(self, particles, binarray):
        """Split particles according to weight.  New particles must have their
        lineage information set correctly and be put into the 
        `particles_created` collection.  The source particle must be added
        to the `particles_split` collection.
        """                
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
                    
                    self.particles_split.add(particle)                    
                    bin.discard(particle)
                    particles.discard(particle)
                    
                    self.particles_created.update(new_particles)
                    bin.update(new_particles)
                    particles.update(new_particles)
    
    def merge_particles(self, particles, binarray):
        """Merge particles according to weight.  The conglomerate particle
        must have its lineage information set properly and be placed into the 
        `particles_created` collection.  The source particles must be added to
        the `particles_merged` collection.
        """
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
                    
                    self.particles_created.add(glom)
                    bin.add(glom)
                    particles.add(glom)
                    
                    self.particles_merged.update(glom_src)
                    self.particles_merged.remove(glom.p_parent)
                    bin.difference_update(glom_src)
                    particles.difference_update(glom_src)
                # end if glom_src
                
                mergable_particles = [particle for particle in bin
                                      if particle.weight <= min_weight]
            # end while len(mergable_particles) > 1
                
    def run_we(self, particles):
        """Run the weighted ensemble bin/split/merge algorithm on the given 
        particles.  Sets the following instance variables (all sets):
          `particles_created`
              particles created for any reason (split or merge)
          `particles_escaped`
              particles which exited the bin space (can be used for rudimentary
              probability recycling)
          `particles_merged`
              particles which ceased to be in existence due to a merge
              (*including* the particle whose phase space point was selected
              as that of the new conglomerate particle)
          `particles_split`
              particles which ceased to be in existence due to a split,
              *including* the parent particle
              
        Returns the set of particles (with properly-set lineage information)
        to be propagated next. Particles in this set without any lineage
        information indicate that their probability is to be recycled.
        """
        assert(particles)
        
        self.particles_created = set()
        self.particles_escaped = set()
        self.particles_merged  = set()
        self.particles_split   = set()

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
            
        # By this point, the self.particles_* collections are properly populated
        # and can (must) be used to assign particle lineage to those particles
        # left untouched by the various WE routines

        log.info('%d particles escaped' % len(self.particles_escaped))
        log.info('%d particles split' % len(self.particles_split))
        log.info('%d particles merged' % len(self.particles_merged))        
        log.info('%d new particles created' % len(self.particles_created))
                 
        # Various sanity checks
        assert abs(particles.norm - 1) < 1.0e-14
        
        for particle in self.particles_merged:
            assert particle not in particles
            
        for particle in self.particles_split:
            assert particle not in particles
        
        for particle in self.particles_created:
            assert particle.p_parent is not None
            assert particle in particles
        
        for particle in self.particles_escaped:
            assert particle in particles
            # Clear particle lineage information so that the particles start
            # new trajectories next round
            particle.parents = []
            particle.p_parent = None
            
        self.continued_particles = set(particles) \
                                 - self.particles_created \
                                 - self.particles_escaped
        
        # Assign lineage information for trajectory continuation
        for particle in self.continued_particles:
            particle.parents = [particle]
            particle.p_parent = particle
        
        # Calculate per-bin population and flux
        last_population = self.bins_population
        last_nparticles = self.bins_nparticles
        self.bins_population = bins.population_array()
        self.bins_nparticles = bins.nparticles_array()
        
        try:
            self.bins_popchange = last_population - self.bins_population
        except (TypeError,ValueError):
            self.bins_popchange = None
            
        self.current_iteration += 1
        
        return particles
