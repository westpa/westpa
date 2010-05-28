from __future__ import division
__metaclass__ = type

import logging
log = logging.getLogger('wemd.we_drivers.we_driver')

import math, random
from copy import copy
import numpy 

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

	self.sim_config = sim_config
	self.initial_pcoord = numpy.array(sim_config.get_list('wemd.initial_pcoord', type=float))
	self.target_pcoord_lower = numpy.array(sim_config.get_list('wemd.target_pcoord_lb', type=float))
	self.target_pcoord_upper = numpy.array(sim_config.get_list('wemd.target_pcoord_ub', type=float))
                            
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

        target_lb = self.target_pcoord_lower
        target_ub = self.target_pcoord_upper

        for particle in particles:

            particle_escaped = False
            assert(len(particle.pcoord) == ndim)
            index = [None] * ndim
            for idim in xrange(0, ndim):
                #index[idim] = numpy.digitize( particle.pcoord[idim], boundaries[idim] )
                for ibound in xrange(0, len(boundaries[idim])-1):
                    if boundaries[idim][ibound] <= particle.pcoord[idim] < boundaries[idim][ibound+1]:
                        index[idim] = ibound

            for idim in xrange(0, ndim):
                if (particle.pcoord[idim] < target_lb[idim]) or (particle.pcoord[idim] > target_ub[idim]):
                    break
                elif idim == (ndim - 1):
                    self.particles_escaped.add(particle)
                    particle_escaped = True

            index = tuple(index)
            if None in index:
                raise ValueError("Error, particle not in binspace %r" % particle.pcoord)

            if particle_escaped == False:
                binarray.bins[index].add(particle)
 
    def split_particles(self, particles, binarray):
        """Split particles according to weight.  New particles must have their
        lineage information set correctly and be put into the 
        `particles_created` collection.  The source particle must be added
        to the `particles_split` collection.
        """             

        for bin in binarray:
            bnorm = bin.norm

            cur_num = len(bin) #current number of particles
            if( (cur_num == 0) or (cur_num >= bin.ideal_num) ): #no need to split
                continue

            particles_list = list(bin)[:]
            weight_key = (lambda p: p.weight)

            #sort particles by weight
            sorted_particles = sorted( particles_list, key = weight_key, reverse = True ) 
            
            proposed_split = [] #contains a list of particles that will be split
                                #[index_into_sorted_particles, nparticle to create, new_weight, orig_weight]

            for i in range(0, len(sorted_particles)):
                proposed_split.append([i, 1, sorted_particles[i].weight, sorted_particles[i].weight ]) #add the particles

            while( cur_num != (bin.ideal_num) ):

                for i in range(0, len(proposed_split)): #loop through proposed_split

                    if( i < (len(proposed_split)-1) ): #if there are more particles, split by next weight
                        while 1: #do { } while ( ... )

                            proposed_split[i][1] += 1 #number of particles to create
                            cur_num += 1
                            proposed_split[i][2] = proposed_split[i][3] / proposed_split[i][1] #update weight

                            if( cur_num == (bin.ideal_num) ):
                                break

                            if( i > 0 ): #more evenly split by checking to see how previous weight compares to this one

                                #weight if we were to split the prev particle again
                                prev_weight = proposed_split[i-1][3] / (proposed_split[i-1][1] + 1)

                                if( prev_weight >= proposed_split[i][2] ): #split particle if it will be greater than the next weight
                                    proposed_split[i-1][1] += 1 #number of particles to create
                                    cur_num += 1
                                    proposed_split[i-1][2] = proposed_split[i-1][3] / proposed_split[i-1][1] #update weight
                                    
                                    if( cur_num == (bin.ideal_num) ):
                                        break
                                
                            if ( proposed_split[i][2] < proposed_split[i+1][2] ): #continue until proposed weight is less than next weight
                                break

                    else: #last particle in the bin split it (for cases with one particle, etc)
                        proposed_split[i][1] += 1 #number of particles to create
                        cur_num += 1
                        proposed_split[i][2] = proposed_split[i][3] / proposed_split[i][1] #update weight

                    if( cur_num == (bin.ideal_num) ):
                        break                    

            x = 0
            for i in range(0, len(proposed_split)):

                particle = sorted_particles[proposed_split[i][0]]
                new_num = proposed_split[i][1]

                if( new_num < 2 ): #no need to do anything to this particle
                    continue

                P = particle.__class__
                assert(new_num >= 2)
                    
                new_weight = particle.weight / new_num
                new_particles = [P(weight = new_weight,
                                   pcoord = copy(particle.pcoord),
                                   p_parent = particle,
                                   parents = [particle])]
    
                new_particles.extend(P(weight = new_weight,
                                       pcoord = copy(particle.pcoord),
                                       p_parent = particle,
                                       parents = [particle])
                                     for i in xrange(1, new_num))
                    
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

            cur_num = len(bin) #current number of particles
            if( (cur_num < 2) or (cur_num <= bin.ideal_num) ): #no need to merge
                continue

            proposed_merge = [] #a list of proposed particles to mnerge,
                                #each element contains an array of particles that will be merged

            mergable_particles = sorted( [particle for particle in bin], key = weight_key )
            for particle in mergable_particles:
                proposed_merge.append([particle])

            while( len(proposed_merge) != bin.ideal_num ):

                #the length of proposed_merge changes
                #so we can't use for i in range(0,len(proposed_merge))
                i = 0
                while i < len(proposed_merge) - 1:

                    if( len(proposed_merge) == bin.ideal_num ):
                        break

                    #len(proposed_merge) changes, i does not
                    while i < len(proposed_merge) - 1:

                        proposed_merge[i].append( proposed_merge[i+1][0] )
                        proposed_merge[i+1].pop(0)

                        if( len(proposed_merge[i+1]) == 0 ):#all elements merged
                            proposed_merge.pop(i+1) #remove empty container

                        new_weight = 0.0
                        for x in range(0, len(proposed_merge[i])):
                            new_weight += proposed_merge[i][x].weight
    
                        if( len(proposed_merge) == bin.ideal_num ):
                            break

                        #weight high, move onto next proposed particle to merge
                        if( new_weight > ideal_weight ):
                            break

                    i += 1
               
            for i in range(0, len(proposed_merge)):

                if len(proposed_merge[i]) > 1: #need 2 particles to merge
                    glom_src = proposed_merge[i]
                else:
                    continue

                glom_weight = 0.0 #get the weight of the new particle
                for x in range(0, len(proposed_merge[i])):
                    glom_weight += proposed_merge[i][x].weight

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
        log.info('%d particles escaped' % len(self.particles_escaped))
        for particle in self.particles_escaped:
            assert particle in particles
            # Clear particle lineage information so that the particles start
            # new trajectories next round
            particle.parents = []
            particle.p_parent = None
            particle.pcoord = copy(self.initial_pcoord)
	self.distribute_particles(particles, bins)

        
        # endpoint bin assignments stored here
        
        # Split/merge particles
        if self.current_iteration > 0:
            self.split_particles(particles, bins)
            self.merge_particles(particles, bins)
            
        # By this point, the self.particles_* collections are properly populated
        # and can (must) be used to assign particle lineage to those particles
        # left untouched by the various WE routines
        log.info('%d particles split' % len(self.particles_split))
        log.info('%d particles merged' % len(self.particles_merged))        
        log.info('%d new particles created' % len(self.particles_created))
        log.info('%d particles total' % len(particles))
        # Various sanity checks
        #assert abs(particles.norm - 1) < 1.0e-14
        
        for particle in self.particles_merged:
            assert particle not in particles
            
        for particle in self.particles_split:
            assert particle not in particles
        
        for particle in self.particles_created:
            assert particle.p_parent is not None
            assert particle in particles
                    
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
