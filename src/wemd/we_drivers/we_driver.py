from __future__ import division
__metaclass__ = type

import logging
log = logging.getLogger('wemd.we_drivers.we_driver')

from itertools import izip
import math, random
from copy import copy
import numpy 

class WEDriver:
    def __init__(self):
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
        
    def sim_init(self, sim_config, sim_config_src):
        self.sim_config = sim_config
        
        sim_config_src.require('wemd.initial_pcoord')
        self.initial_pcoord = sim_config['wemd.initial_pcoord'] \
                            = numpy.array(sim_config_src.get_list('wemd.initial_pcoord', type=float))

    def runtime_init(self, runtime_config):
        self.runtime_config = runtime_config
        
    def make_bins(self):
        """Create an array of ParticleCollection objects appropriate to this WE
        method."""
        raise NotImplementedError
                                        
    def distribute_particles(self, particles, binarray):
        """Assign particles to bins.  Particles cannot fall out of the bin space.
        """
        
        from numpy import digitize
         
        boundaries = binarray.boundaries
        ndim = binarray.ndim

        n_targets = len(self.target_pcoords)
        target_bounds = numpy.empty((n_targets, ndim, 4), numpy.float64)
        for itarg in xrange(0,n_targets):
            target_bounds[itarg,:,0] = float('-inf')
            target_bounds[itarg,:,1] = self.target_pcoords[itarg,:,0]
            target_bounds[itarg,:,2] = self.target_pcoords[itarg,:,1]
            target_bounds[itarg,:,3] = float('inf')
              
        particle_list = numpy.empty((len(particles),), numpy.object_)
        recycled_particles = set()
        particle_list[:] = list(particles)
        #print particle_list, type(particle_list)
        pcoords = numpy.empty((len(particles), ndim), numpy.float64)
        for (ipart, particle) in enumerate(particle_list):
            pcoords[ipart] = particle.pcoord
        indices = numpy.empty((pcoords.shape[0], ndim), numpy.uintp)
        
        # Find particles at the sink(s)
        for itarget in xrange(0, target_bounds.shape[0]):
            for idim in xrange(0, ndim):
                #print target_bounds[itarget,idim,:]
                indices[:, idim] = numpy.digitize(pcoords[:, idim], target_bounds[itarget,idim,:])
            # Every row which is all 2s indicates a particle which is entirely within the given target region
            
            recycled_particles.update(particle_list[numpy.all(indices==2, axis=-1)])
        
        # Assign particles to bins
        particles = particles - recycled_particles            
        particle_list = numpy.empty((len(particles),), numpy.object_)
        particle_list[:] = list(particles)
        pcoords = numpy.empty((len(particles), ndim), numpy.float64)
        for (ipart, particle) in enumerate(particle_list):
            pcoords[ipart] = particle.pcoord
        indices = numpy.empty((pcoords.shape[0], ndim), numpy.uintp)    
        
        for idim in xrange(0, ndim):
            indices[:, idim] = numpy.digitize(pcoords[:, idim], boundaries[idim])
            if (indices[indices[:, idim] == 0].size
               or indices[indices[:, idim] == len(boundaries[idim])].size):
                raise ValueError('particle not in bin space')
        # Rebase so that we have a valid index into the array of particle collections
        indices -= 1

        for (particle, index) in izip(particle_list, indices):
            #print particle, index
            binarray.bins[tuple(index)].add(particle)
        self.particles_escaped.update(recycled_particles)
        
    def split_particles_weight(self, particles, binarray):
        """Split particles according to weight.  New particles must have their
        lineage information set correctly and be put into the 
        `particles_created` collection.  The source particle must be added
        to the `particles_split` collection.
        """                
        for bin in binarray:
            bnorm = bin.norm
            ideal_weight = bnorm / bin.ideal_num
            for particle in list(bin):

                if particle.weight > 2.0*ideal_weight:
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
    
    def merge_particles_weight(self, particles, binarray):
        """Merge particles according to weight.  The conglomerate particle
        must have its lineage information set properly and be placed into the 
        `particles_created` collection.  The source particles must be added to
        the `particles_merged` collection.
        """
        weight_key = (lambda p: p.weight)
        
        for bin in binarray:
            ideal_weight = bin.norm / bin.ideal_num
            min_weight = ideal_weight * 0.5
            max_weight = ideal_weight * 1.5
            
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
                    
    def split_particles_number(self, particles, binarray):
        """Split particles according to weight, while making the number of particles
        in each bin eqaul to the ideal number in the bin.
        New particles must have their lineage information set correctly and be put 
        into the `particles_created` collection.  The source particle must be added
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


    def merge_particles_number(self, particles, binarray):
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

    def split_merge_adjust(self):
        """Fixes up parents of particles that were split then merged (etc) to prevent
        phantom particles from showing up.
        """
        particles_created_discard = []
        particles_merged_discard = []
        particles_split_discard = []
        for particle in self.particles_created: 
            pp_remove = []
            pp_append = []

            for i in xrange(0, len(particle.parents)):
                pp = particle.parents[i]
                if pp in self.particles_created:
                    particles_created_discard.append(pp) #temporary particle not in db
                    
                    if pp in self.particles_merged:
                        pp_append.append(pp.p_parent)
                        pp_remove.append(pp)
                        particles_merged_discard.append(pp)
                        if( pp.particle_id == particle.p_parent.particle_id ): 
                            #print("replacing p_parent merged")
                            particle.p_parent = pp.p_parent
                            
                        assert pp.p_parent not in self.particles_created
         
                    elif pp in self.particles_split:                
                        pp_append.append(pp.p_parent)
                        pp_remove.append(pp)
                        particles_split_discard.append(pp)
                        if( pp.particle_id == particle.p_parent.particle_id ):
                            particle.p_parent = pp.p_parent

                        assert pp.p_parent not in self.particles_created
                    else: #self.particles_merged.remove(glom.p_parent) <-- so it was merged but is not in particles_merged                        
                        pp_append.append(pp.p_parent)
                        pp_remove.append(pp)
                        particles_merged_discard.append(pp)
                        if( pp.particle_id == particle.p_parent.particle_id ): 
                            particle.p_parent = pp.p_parent
                            
                        assert pp.p_parent not in self.particles_created
                    
            for pp in pp_remove:
                particle.parents.remove(pp)   
            for pp in pp_append:
                particle.parents.append(pp)              

        for p in particles_created_discard:
            self.particles_created.discard(p)
        for p in particles_merged_discard:
            self.particles_merged.discard(p)
        for p in particles_split_discard:
            self.particles_split.discard(p)   
              
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
        ndim = bins.ndim
        source_pcoords = self.source_pcoords

        #use one coordinate recycling
        if not source_pcoords:
            # Clear particle lineage information so that the particles start
            # new trajectories next round
            for particle in self.particles_escaped:
                assert particle in particles                
                particle.parents = []
                particle.p_parent = None
                particle.pcoord = copy(self.initial_pcoord)
        #recycle proportional to the weight in a specified region
        #source_pcoords['Region Name'] is a dictionary with the keys 'pcoord', 'region', and 'weight'
        else:
            npcoords = len(source_pcoords)
            kcoord = source_pcoords.keys()
                        
            for particle in self.particles_escaped:
                assert particle in particles
                
                region_prob = numpy.zeros((npcoords,),dtype=numpy.float64)
                #use initial weight
                for icoord in xrange(0,npcoords):      
                    region_prob[icoord] = source_pcoords[kcoord[icoord]]['weight']

                icoord = 0
                while(True):

                    if numpy.random.random_sample() < region_prob[icoord]:                       
                        particle.parents = []
                        particle.p_parent = None
                        particle.pcoord = copy(source_pcoords[kcoord[icoord]]['pcoord'])
                        particle.initial_region = kcoord[icoord]
                        break
                                    
                    icoord = (icoord+1) % npcoords

        self.distribute_particles(particles, bins)
        # endpoint bin assignments stored here
        
        # Split/merge particles
        if self.current_iteration > 0:
            self.split_particles_weight(particles, bins)
            self.merge_particles_weight(particles, bins)
            self.split_merge_adjust()
            self.split_particles_number(particles, bins)
            self.split_merge_adjust()
            self.merge_particles_number(particles, bins)
            self.split_merge_adjust()

          
        for bin in bins:
            ideal_weight = bin.norm / bin.ideal_num
            for particle in bin:
                if particle.weight < (0.1 * ideal_weight) or (particle.weight > 4.0 * ideal_weight):
                    log.warning('Warning: particle weight out of range:\n\
                                 [%r , %r]: %r' % (ideal_weight * 0.1, ideal_weight * 4.0, particle.weight))
                     
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
