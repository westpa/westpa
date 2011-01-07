import logging
log = logging.getLogger('wemd.we_drivers.fixed_bins_particles')

import numpy
import math, random
from copy import copy

from we_driver import WEDriver
from fixed_bins import FixedBinWEDriver
from wemd.core.binarrays import Bin, BinArray
from wemd.core import ConfigError

class FixedBinParticlesWEDriver(FixedBinWEDriver):
    def __init__(self, sim_config):
        super(FixedBinParticlesWEDriver,self).__init__(sim_config)
        self.configure_bin_params(sim_config)
        
        # List to keep track of all particles that were created in split and merge
        self.particles_created_all = set()
    
    def split_particles(self, particles, binarray):
        """ This routine splits particles in a bin until the number of particles in the bin is
        greater than or equal to the the ideal number of particles. Particles are split in descending 
        order ranked by their weight.
        """
        weight_key = (lambda p: p.weight)
        
        for bin in binarray:
            bnorm = bin.norm
            ideal_weight = bnorm / bin.ideal_num
            ntot = len(bin)
            bin_temp = list(bin)[:]
            bin_ordered = sorted(bin_temp,key = weight_key,reverse = True)
            if len(bin) != 0 and len(bin) < bin.ideal_num: 
                for particle in bin_ordered:
                    m = int(math.ceil(particle.weight / ideal_weight))
                
                    P = particle.__class__
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
                    
                    self.particles_created_all.update(new_particles)
                
                    ntot += (m - 1)
                
                    if ntot >= bin.ideal_num:
                        break
                        
            
            if(len(bin)):
                assert(len(bin) >= bin.ideal_num), "Split resulted in less particles than ideal number"        
                
    def merge_particles(self, particles, binarray):
        """ This routine merges the two particles in the bin with the smallest weight iteratively, until 
        the number of particles in the bin is equal to the ideal number of particles
        """
        weight_key = (lambda p: p.weight)

        for bin in binarray:
            if len(bin) > 0:
                assert(len(bin) >= bin.ideal_num), "Less than the ideal number of particles preceeding merge"
            ideal_weight = bin.norm / bin.ideal_num
            
            while(len(bin) > bin.ideal_num):
                a = sorted(list(bin), key = weight_key)[:2]
                glom_src = []
                glom_weight = 0.0
                for particle in a:
                    glom_src.append(particle)
                    glom_weight += particle.weight
                
                P = glom_src[0].__class__
                glom = P()
                u = random.uniform(0,glom_weight)
                
                if(u < glom_src[0].weight):
                    glom.pcoord = glom_src[0].pcoord
                    if glom_src[0] in self.particles_created_all:
                        glom.p_parent = glom_src[0].p_parent
                    else:
                        glom.p_parent = glom_src[0]
                else:
                    glom.pcoord = glom_src[1].pcoord
                    if glom_src[1] in self.particles_created_all:
                        glom.p_parent = glom_src[1].p_parent
                    else:
                        glom.p_parent = glom_src[1]

                    
                glom.weight = glom_weight
                
                # Properly handle parents if glom_src particles were created during split/merge
                if glom_src[0] in self.particles_created_all:
                    glom.parents.update(glom_src[0].parents)
                else:
                    glom.parents.add(glom_src[0])
                   
                if glom_src[1] in self.particles_created_all:
                    glom.parents.update(glom_src[1].parents)
                else:
                    glom.parents.add(glom_src[1])
                                    
                glom.particle_id = glom.p_parent.particle_id

                self.particles_created.add(glom)
                bin.add(glom)
                particles.add(glom)

                self.particles_created_all.add(glom)

                self.particles_merged.update(glom_src)
                bin.difference_update(glom_src)
                particles.difference_update(glom_src)
                
                # Remove glom_src particles from particles_created list if they are merged into glom
                self.particles_created.difference_update(glom_src)
                    
            log.debug("Post-merge Bin: %d Norm: %g Length: %d" % (bin.bin_id,bin.norm,len(bin)))
        self.particles_created_all = set()
            
                
                
    