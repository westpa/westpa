from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)
import numpy
import operator
from itertools import izip, chain
from math import ceil
from copy import copy
import random

from collections import namedtuple
from wemd.util.miscfn import vgetattr
from wemd import Segment

RecyclingInfo = namedtuple('RecyclingInfo', ['count', 'weight'])

class WEMDWEDriver:
    weight_split_threshold = 2.0
    weight_merge_cutoff = 1.0
        
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
        # (count recycled, probability recycled) tuples
        self.recycle_from = list()
        
        # (count_recycled, probability recycled) lists
        self.recycle_to = list()
                        
    def run_we(self, segments, region_set):
        '''Run the weighted ensemble algorithm on the given segments, binning on their final
        positions.  The endpoint_type field of each segment is updated appropriately.  Returns the
        new set of segments to be propagated.'''
        
        segments = list(segments)
        
        # Everything continues unless it gets split, merged, or recycled
        for segment in segments:
            segment.endpoint_type = Segment.SEG_ENDPOINT_TYPE_CONTINUES

        self.recycle_from = [[0, 0.0] for istate in xrange(0, len(self.sim_manager.system.target_states))]
        self.recycle_to   = [[0, 0.0] for istate in xrange(0, len(self.sim_manager.system.initial_states))]  
                        
        # Distribute segments into bins based on their endpoints 
        region_set.clear()
        endpoint_coords = [segment.pcoord[-1] for segment in segments]
        for (segment, bin) in izip(segments, region_set.map_to_bins(endpoint_coords)):
            bin.add(segment)
            
        bins = numpy.array(region_set.get_all_bins(), numpy.object_)
        
        # Determine which segments terminate in recycling events
        self.recycle_particles(segments, self.sim_manager.system.target_states)

        # Sanity check
        bin_counts = vgetattr('count', bins)
        target_counts = vgetattr('target_count', bins)            
        if bin_counts[target_counts == 0].any():
            raise ValueError('simulation in forbidden region')
        
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
                
        # Replace these simple lists with named tuples to make others' lives easier        
        self.recycle_to = [RecyclingInfo(*dest) for dest in self.recycle_to]
        self.recycle_from = [RecyclingInfo(*src) for src in self.recycle_from]
            
        new_segments = []
        new_pcoord_array = self.sim_manager.system.new_pcoord_array
        for segment in chain(*bins):
            if log.isEnabledFor(logging.DEBUG):
                log.debug('examining segment {!r}'.format(segment))
            if segment.seg_id is not None:
                # A simple continuation; we need to create a new segment for it here
                new_segment = Segment(weight = segment.weight,
                                      pcoord = new_pcoord_array(),
                                      p_parent_id = segment.seg_id,
                                      parent_ids = [segment.seg_id])
                new_segment.pcoord[0] = segment.pcoord[-1]
                new_segments.append(new_segment)
            else:
                # A segment resulting from recycling, splitting, or merging;
                # There is no need to create a new segment
                assert segment.p_parent_id is not None
                assert segment.p_parent_id in segment.parent_ids
                assert segment.pcoord is not None
                new_segments.append(segment)
        return new_segments
    
    def _log_bin_stats(self, bin, heading=None, level=logging.DEBUG):
        if log.isEnabledFor(level):
            weights = numpy.array(map(operator.attrgetter('weight'), bin))
            weights.sort()
            bin_label = getattr(bin, 'label', None) or ''
            log_fmt = '\n      '.join(['', 
                                         'stats for bin {bin_label!r} {heading}',
                                         '  count: {bin.count:d}, target count: {bin.target_count:d}',
                                         '  total weight: {bin.weight:{weight_spec}},  ideal weight: {ideal_weight:{weight_spec}}',
                                         '  mean weight:  {mean_weight:{weight_spec}},  stdev weight: {stdev_weight:{weight_spec}}',
                                         '  min weight:   {min_weight:{weight_spec}},  med weight  : {median_weight:{weight_spec}}'
                                         +', max weight: {max_weight:{weight_spec}}'])
            log_msg = log_fmt.format(log_fmt, weight_spec='<12.6e', 
                                     bin_label=bin_label, heading=heading, bin=bin,
                                     ideal_weight = bin.weight/bin.target_count,
                                     mean_weight = weights.mean(), stdev_weight = weights.std(),
                                     min_weight = weights[0], median_weight=numpy.median(weights), max_weight=weights[-1])
            log.log(level, log_msg)
        
    def recycle_particles(self, segments, target_states):
        '''Create new particles with the weight of the given particles whose progress coordinates are 
        chosen probabilistically from the set of initial states.'''
        
        system = self.sim_manager.system
        for (itarget, target_state) in enumerate(target_states):
            bin = target_state.bin
            
            if log.isEnabledFor(logging.DEBUG):
                log.debug('{0.count:d} particles recycled from bin {0!r}'.format(bin))
            
            if bin.count == 0:
                continue
            
            self.recycle_from[itarget] = (bin.count, bin.weight)
            
            new_segments = []
            for segment in bin:
                segment.endpoint_type = Segment.SEG_ENDPOINT_TYPE_RECYCLED
                new_segment = Segment(weight = segment.weight,
                                      pcoord = system.new_pcoord_array())
                new_segments.append(new_segment)
            
            rstates = [(istate, state) for (istate,state) in enumerate(system.initial_states) if state.recycle_prob > 0]
            
            # Sequentially add probabilities for initial states.
            state_disc = numpy.add.accumulate([state.recycle_prob for (istate,state) in rstates])
            
            # Generate one random number for each segment being recycled from this bin
            # Bin on (0, state_disc[-1]) in case probabilities don't sum to exactly one
            irstates = numpy.digitize(numpy.random.uniform(0, state_disc[-1], size=bin.count),state_disc)
            
            for (new_segment,irstate) in izip(new_segments,irstates):
                iistate, istate = rstates[irstate]
                p_parent_id = -(iistate+1)
                new_segment.pcoord[0] = istate.pcoord[:]
                new_segment.p_parent_id = p_parent_id
                new_segment.parent_ids = set((p_parent_id,))
                self.recycle_to[iistate][0] += 1
                self.recycle_to[iistate][1] += new_segment.weight
                istate.bin.add(new_segment)
            bin.clear()
            
    def split_by_weight(self, bin):
        '''Split overweight particles'''
        target_count = bin.target_count
        ideal_weight = bin.weight / target_count
        
        segments = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
        weights = numpy.array(map(operator.attrgetter('weight'), segments))
        
        for segment in segments[weights > self.weight_split_threshold * ideal_weight]:
            segment.endpoint_type = Segment.SEG_ENDPOINT_TYPE_CONTINUES
            bin.remove(segment)
            m = int(ceil(segment.weight / ideal_weight))
            new_segments = []
            new_weight = segment.weight / m
            if segment.seg_id is not None:
                new_pcoord = segment.pcoord[-1]
            else:
                new_pcoord = segment.pcoord[0]
            new_parent_id = segment.seg_id if segment.seg_id is not None else segment.p_parent_id
            for inew in xrange(0,m):
                new_segment = Segment(weight = new_weight,
                                      p_parent_id = new_parent_id,
                                      parent_ids = [new_parent_id],
                                      pcoord = self.sim_manager.system.new_pcoord_array())
                new_segment.pcoord[0] = new_pcoord
                new_segments.append(new_segment)
            
            if log.isEnabledFor(logging.DEBUG):
                log.debug('splitting {!r} into {:d}: {!r}'.format(segment, len(new_segments), new_segments)) 
            bin.update(new_segments)
    
    def merge_by_weight(self, bin):
        '''Merge underweight particles'''
        target_count = bin.target_count
        ideal_weight = bin.weight / target_count
        
        segments = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
        weights = numpy.array(map(operator.attrgetter('weight'), segments))
        cumul_weight = numpy.add.accumulate(weights)
        
        to_merge = segments[cumul_weight <= ideal_weight*self.weight_merge_cutoff]
        if len(to_merge) < 2:
            return
        else:
            to_merge_ids = [segment.seg_id if segment.seg_id is not None else segment.p_parent_id for segment in to_merge]
        
        assert cumul_weight[len(to_merge)-1] == sum(weights[:len(to_merge)])
        glom = Segment(weight = cumul_weight[len(to_merge)-1],
                       parent_ids = set(to_merge_ids),
                       pcoord = self.sim_manager.system.new_pcoord_array()) 
        iparent = numpy.digitize((random.uniform(0,glom.weight),),cumul_weight[:len(to_merge)])[0]
        if to_merge[iparent].seg_id is not None:
            glom.pcoord[0] = to_merge[iparent].pcoord[-1]
        else:
            glom.pcoord[0] = to_merge[iparent].pcoord[0]
        glom.p_parent_id = to_merge_ids[iparent]
        
        if log.isEnabledFor(logging.DEBUG):
            log.debug('merging {:d} {!r} into {!r}'.format(len(to_merge), to_merge, glom))
        
        # The primary parent is marked as continuing; others are marked as ending in a merge
        for (iseg, segment) in enumerate(to_merge):
            if iseg != iparent:
                segment.endpoint_type = Segment.SEG_ENDPOINT_TYPE_MERGED
            else:
                segment.endpoint_type = Segment.SEG_ENDPOINT_TYPE_CONTINUES
                
        bin.difference_update(to_merge)
        bin.add(glom)
    
    def adjust_count(self, bin):
        '''Adjust the particle count of a bin (by splits/merges) so that it exactly equals the target count'''
        weight_getter = operator.attrgetter('weight')
        
        # split
        while bin.count < bin.target_count:
            # Always split the highest probability particle into two
            particles = list(sorted(bin, key=weight_getter))
            parent = particles[-1]
            # If the highest-weight particle is newly-created, it will not have a seg_id;
            # In that case, we essentially replace that segment with two new ones
            children = [None] * 2
            for i in xrange(0,2):
                child = Segment(weight = parent.weight/2,
                                pcoord = self.sim_manager.system.new_pcoord_array())
                if parent.seg_id is not None:
                    # This is a simple continuation
                    child.p_parent_id = parent.seg_id
                    child.parent_ids = set([parent.seg_id])
                    child.pcoord[0] = parent.pcoord[-1]
                else:
                    # The parent segment is the result of a recycling, split, or merge
                    child.p_parent_id = parent.p_parent_id
                    child.parent_ids = set([parent.p_parent_id])
                    child.pcoord[0] = parent.pcoord[0]
                children[i] = child              

            if log.isEnabledFor(logging.DEBUG):
                log.debug('splitting {!r} into {:d}: {!r}'.format(parent, len(children), children))             
            bin.remove(parent)
            bin.update(children)
        
        # merge
        while bin.count > bin.target_count:
            # Always merge the two lowest-probability particles
            segments = list(sorted(bin, key=weight_getter))
            parents = segments[0:2]
            glom_weight = parents[0].weight + parents[1].weight
            parent_ids = [segment.seg_id if segment.seg_id is not None else segment.p_parent_id
                          for segment in parents]
            assert None not in parent_ids
            glom = Segment(weight = glom_weight,
                           parent_ids = set(parent_ids),
                           pcoord = self.sim_manager.system.new_pcoord_array()
                          )
            if random.uniform(0, glom_weight) < parents[0].weight:
                iparent = 0
            else:
                iparent = 1
            
            glom.p_parent_id = parent_ids[iparent]
            iother = 1 if iparent == 0 else 0
            parents[iparent].endpoint_type = Segment.SEG_ENDPOINT_TYPE_CONTINUES
            parents[iother].endpoint_type  = Segment.SEG_ENDPOINT_TYPE_MERGED

            if parents[iparent].seg_id is not None:
                glom.pcoord[0] = parents[iparent].pcoord[-1]
            else:
                glom.pcoord[0] = parents[iparent].pcoord[0]

            if log.isEnabledFor(logging.DEBUG):
                log.debug('merging {:d} {!r} into {!r}'.format(len(parents), parents, glom))
                
            bin.difference_update(parents)
            bin.add(glom)
            
