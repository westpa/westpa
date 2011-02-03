from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)
import numpy
import operator
from itertools import izip, chain
from math import floor, ceil
from copy import copy
import random

from wemd.types import Particle
from wemd.systems import WEMDSystem

class WEMDWEDriver:
    weight_split_threshold = 2.0
    weight_merge_cutoff = 1.0
    
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
        # IDs of segments which terminate in recycling
        self.recycle_terminations = set()
        
        # IDs of segments which are merged into others
        self.merge_terminations = set()
        
    def dump_state(self):
        pass
    
    def restore_state(self, stateinfo):
        pass
        
    def run_we(self, particles, region_set):
        '''Run the weighted ensemble algorithm on the given particles, setting the particle
        collections appropriately'''
        particles = list(particles)
        self.recycle_terminations = set()
        self.merge_terminations = set()

        
        # Unless otherwise manhandled, all particles will merely continue propagation and therefore are 
        # in a 1-to-1 correspondence with the segments from which they were created
        for particle in particles:
            particle.genesis = Particle.GEN_CONTINUATION
                
        # Distribute particles 
        region_set.clear()
        for (particle, bin) in izip(particles, region_set.map_to_bins(particle.pcoord for particle in particles)):
            bin.add(particle)
            
        bins = numpy.array(region_set.get_all_bins(), numpy.object_)
        target_counts = numpy.array(map(operator.attrgetter('target_count'), bins))
                
        # Recycle particles which are in bins with target_count==0 (indicating a sink)
        for bin in bins[target_counts==0]:
            self.recycle_particles(bin)
        
        # Regardless of current particle count, always split overweight particles and merge underweight particles
        # Then and only then adjust for correct particle count
        for (ibin,bin) in enumerate(bins):
            if bin.count == 0: continue
            if bin.label is None: bin.label = str(ibin)
            self._log_bin_stats(bin, 'entering WE')
            self.split_by_weight(bin)
            self._log_bin_stats(bin, 'after split by weight')
            self.merge_by_weight(bin)
            self._log_bin_stats(bin, 'after merge by weight')
            self.adjust_count(bin)
            self._log_bin_stats(bin, 'after count adjustment')
        
        # Make sure no invalid seg_ids slipped through
        assert None not in self.recycle_terminations
        assert None not in self.merge_terminations
        
        # merge_terminations may include negative seg_ids if recycled particles are merged
        # eliminate these
        self.merge_terminations = set(seg_id for seg_id in self.merge_terminations if seg_id >= 0)
            
        # Return all particles, existing or created
        return set(chain(*bins))
    
    def _log_bin_stats(self, bin, heading=None, level=logging.DEBUG):
        if log.isEnabledFor(level):
            weights = numpy.array(map(operator.attrgetter('weight'), bin))
            weights.sort()
            bin_label = getattr(bin, 'label', None)
            bin_label = bin_label or ''
            log.log(level, "stats for bin %r %s" % (bin_label, heading))
            log.log(level, 'count: %d, target count: %d' % (bin.count, bin.target_count))
            log.log(level, 'total weight: %g' % bin.weight)
            log.log(level, 'ideal weight: %g' % (bin.weight / bin.target_count))
            log.log(level, 'mean weight: %g, stdev weight: %g' % (weights.mean(), weights.std()))
            log.log(level, 'min weight: %g, med weight: % g, max weight: %g'
                           % (weights[0], numpy.median(weights), weights[-1]))
            
        
    def recycle_particles(self, bin):
        '''Create new particles with the weight of the given particles whose progress coordinates are 
        chosen probabilistically from the set of initial states.'''
        
        if bin.count == 0:
            log.debug('no particles recycled')
            return
        else:
            log.debug('%d particles recycled from bin %r' % (bin.count,bin))
            
        particles = list(bin)
                        
        system = self.sim_manager.system
        self.recycle_terminations.update(particle.seg_id for particle in particles)
        
        new_particles = [Particle(weight = particle.weight, 
                                  pcoord = numpy.empty((system.pcoord_ndim,), system.pcoord_dtype),
                                  genesis = Particle.GEN_RECYCLE)
                         for particle in particles]
        
        if len(system.initial_distribution) == 1:
            init_pcoord = system.initial_distribution[0][WEMDSystem.INITDIST_PCOORD]
            for new_particle in new_particles:
                new_particle.pcoord = init_pcoord 
                new_particle.source_id = 0
                new_particle.seg_id = -1
            system.initial_distribution[0][WEMDSystem.INITDIST_BIN].update(new_particles)
        else:   
            # Sequentially add probabilities for initial states.
            state_disc = numpy.add.accumulate([state[WEMDSystem.INITDIST_PROB] for state in system.initial_distribution])
            log.debug('state_disc = %r' % state_disc)
            
            # Bin on (0, state_disc[-1]) in case probabilities don't sum to exactly one
            istates = numpy.digitize(numpy.random.uniform(0, state_disc[-1], size=len(particles)),state_disc)
            
            for (new_particle,istate) in izip(new_particles,istates):
                new_particle.pcoord = system.initial_distribution[istate][WEMDSystem.INITDIST_PCOORD]
                new_particle.source_id = istate
                new_particle.seg_id = -(istate+1)
                system.initial_distribution[istate][WEMDSystem.INITDIST_BIN].add(new_particle)
        bin.clear()
            
    def split_by_weight(self, bin):
        '''Split overweight particles'''
        target_count = bin.target_count
        ideal_weight = bin.weight / target_count
        
        particles = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
        weights = numpy.array(map(operator.attrgetter('weight'), particles))
        
        for particle in particles[weights > self.weight_split_threshold * ideal_weight]:
            bin.remove(particle)
            m = int(ceil(particle.weight / ideal_weight))
            bin.update(Particle(genesis = Particle.GEN_SPLIT,
                                weight = particle.weight / m,
                                pcoord = copy(particle.pcoord),
                                p_parent_id = particle.seg_id,
                                parent_ids = set((particle.seg_id,))) for i in xrange(0, m))
    
    def merge_by_weight(self, bin):
        '''Merge underweight particles'''
        target_count = bin.target_count
        ideal_weight = bin.weight / target_count
        
        particles = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
        weights = numpy.array(map(operator.attrgetter('weight'), particles))
        cumul_weight = numpy.add.accumulate(weights)
        
        to_merge = particles[cumul_weight <= ideal_weight*self.weight_merge_cutoff]
        to_merge_ids = [particle.p_parent_id if particle.p_parent_id is not None else particle.seg_id
               for particle in to_merge]
        if len(to_merge) < 2:
            return
        
        assert cumul_weight[len(to_merge)-1] == sum(weights[:len(to_merge)]) 
        glom = Particle(genesis = Particle.GEN_MERGE,
                        weight = cumul_weight[len(to_merge)-1],
                        parent_ids = set(to_merge_ids))
        
        iparent = numpy.digitize((random.uniform(0,glom.weight),),cumul_weight[:len(to_merge)])[0]
        glom.pcoord = copy(to_merge[iparent].pcoord)
        glom.p_parent_id = to_merge_ids[iparent]
        
        bin.difference_update(to_merge)
        self.merge_terminations.update(to_merge_ids[i] for i in xrange(0,len(to_merge)) if i != iparent)
        
        bin.add(glom)
    
    def adjust_count(self, bin):
        '''Adjust the particle count of a bin (by splits/merges) so that it exactly equals the target count'''
        weight_getter = operator.attrgetter('weight')
        
        # split
        while bin.count < bin.target_count:
            # Always split the highest probability particle into two
            particles = list(sorted(bin, key=weight_getter))
            parent = particles[-1]
            # If we are splitting a particle that already has been processed, inherit its parent ID as
            # the parent; otherwise, set the parent ID as normal
            child_pid = parent.p_parent_id if parent.p_parent_id is not None else parent.seg_id
            children = [Particle(genesis = Particle.GEN_SPLIT,
                                 weight = parent.weight / 2,
                                 pcoord = copy(parent.pcoord),
                                 p_parent_id = child_pid,
                                 parent_ids = set((child_pid,))) for i in xrange(0,2)]
            bin.remove(parent)
            bin.update(children)
        
        # merge
        while bin.count > bin.target_count:
            # Always merge the two lowest-probability particles
            particles = list(sorted(bin, key=weight_getter))
            parents = particles[0:2]
            glom_weight = parents[0].weight + parents[1].weight
            parent_ids = [p.p_parent_id if p.p_parent_id is not None else p.seg_id
                          for p in parents]
            assert None not in parent_ids
            glom = Particle(genesis = Particle.GEN_MERGE,
                            weight = glom_weight,
                            parent_ids = set(parent_ids),
                            )
            if random.uniform(0, glom_weight) < parents[0].weight:
                glom.p_parent_id = parent_ids[0]
                glom.pcoord = copy(parents[0].pcoord)
                self.merge_terminations.add(parent_ids[1])
            else:
                glom.p_parent_id = parent_ids[1]
                glom.pcoord = copy(parents[1].pcoord)
                self.merge_terminations.add(parent_ids[0])
                
            bin.difference_update(parents)
            bin.add(glom)
            
