from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)
import numpy
import operator
from math import ceil
import random

import wemd
from wemd import Segment

class WEMDWEDriver:
    weight_split_threshold = 2.0
    weight_merge_cutoff = 1.0
        
    def __init__(self):
        self.system = wemd.rc.get_system_driver()
        self.do_adjust_counts = wemd.rc.config.get_bool('we.adjust_counts', True)
        log.info('Adjust counts to exactly match target_counts: {}'.format(self.do_adjust_counts))

        # RegionSet containing binned segments from concluded iteration (reset on entry to run_we)
        # These segments will have their endpoint_type field set appropriately on exit from run_we()
        self.completed_segs_rset = None
        
        # A mapping of seg_id to the corresponding segment among those from self.completed_segs_rset
        self.completed_segs_map = None
        
        # RegionSet containing binned segments for new iteration, with parent_ids (corresponding to segments
        # in self.completed_segs_rset set) describing how the new segments were constructed (i.e. continuation/split/merge)
        # These segments have their initpoint_type field set appropriately on exit from run_we()
        # reset on entry to run_we()
        self.new_iter_rset = None
        
        
                        
    def run_we(self, region_set, new_trajectory_segs):
        '''Run weighted ensemble split/merge on the segments in the given ``region_set``, setting their
        ``endpoint_type`` field appropriately.  The region set must not contain any segments in regions
        with zero target count, so recycling must be handled prior to calling this function. To 
        assist in recycling probability to newly-created trajectories, the given ``new_trajectory_segs``
        are assumed to start new trajectories for the next iteration, and will be considered in particles-per-bin
        and weight calculations. 
        
        Note that clever use of ``new_trajectory_segs`` could be used to inject new
        trajectories into a running weighted ensemble simulation -- e.g. simulating pulses of reactants being
        added -- if proper normalization (sum of weights in region_set + sum of weights in new_trajectory_segs = 1)
        occurs prior to calling this function.
        '''
        
        del self.completed_segs_rset, self.new_iter_rset, self.completed_segs_map
        self.completed_segs_rset = region_set
        self.new_iter_rset = self.system.new_region_set()
        completed_segs_map = self.completed_segs_map = dict()
        
        prev_iter_bins = region_set.get_all_bins()
        next_iter_bins = self.new_iter_rset.get_all_bins()
        
        # Check to make sure that no simulations exist where there should not be any.
        for bin in region_set.bins:
            if bin.target_count == 0 and bin.weight > 0:
                raise ValueError('a forbidden region ({!r}) is populated'.format(bin))
        
        # Assume until splits/merges that everything from the last propagation phase will continue in the next
        # propagation phase.
        for ibin, bin in enumerate(prev_iter_bins):
            for segment in bin:
                segment.endpoint_type = Segment.SEG_ENDPOINT_CONTINUES
                
                new_segment = Segment(parent_id=segment.seg_id,
                                      weight=segment.weight,
                                      wtg_parent_ids=[segment.seg_id],
                                      pcoord=self.system.new_pcoord_array())
                new_segment.pcoord[0,:] = segment.pcoord[-1,:]
                
                next_iter_bins[ibin].add(new_segment)
                completed_segs_map[segment.seg_id] = segment
                
        # New trajectory segments should be ready to go, with appropriate values set
        # for initpoint_type
        if __debug__:
            for segment in new_trajectory_segs:
                assert segment.initpoint_type == Segment.SEG_INITPOINT_NEWTRAJ
        if new_trajectory_segs:
            self.new_iter_rset.assign_to_bins(new_trajectory_segs, key=Segment.initial_pcoord)
        
        # Regardless of current particle count, always split overweight particles and merge underweight particles
        # Then and only then adjust for correct particle count
        for (ibin,bin) in enumerate(next_iter_bins):
            if bin.count == 0: continue
            if bin.label is None: bin.label = str(ibin)
            self._log_bin_stats(bin, 'entering WE')
            self.split_by_weight(bin)
            self._log_bin_stats(bin, 'after split by weight')
            self.merge_by_weight(bin)
            self._log_bin_stats(bin, 'after merge by weight')
            if self.do_adjust_counts:
                self.adjust_count(bin)
                self._log_bin_stats(bin, 'after count adjustment')
            
        return self.new_iter_rset
                            
    
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
                    
    def split_by_weight(self, bin):
        '''Split overweight particles'''
        target_count = bin.target_count
        ideal_weight = bin.weight / target_count
        
        segments = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
        weights = numpy.array(map(operator.attrgetter('weight'), segments))
        to_split = segments[weights > self.weight_split_threshold*ideal_weight]
        
        for segment in to_split:
            bin.remove(segment)
            m = int(ceil(segment.weight / ideal_weight))
            new_segments = []
            for inew in xrange(0,m):
                new_segment = Segment(weight = segment.weight / m,
                                      parent_id = segment.parent_id,
                                      wtg_parent_ids = set(segment.wtg_parent_ids),
                                      pcoord = segment.pcoord.copy())
                new_segment.pcoord[0,:] = segment.pcoord[0,:]
                new_segments.append(new_segment)
            
            if log.isEnabledFor(logging.DEBUG):
                log.debug('splitting {!r} into {:d}: {!r}'.format(segment, len(new_segments), new_segments)) 
            bin.update(new_segments)
    
    def merge_by_weight(self, bin):
        '''Merge underweight particles'''
        target_count = bin.target_count
        ideal_weight = bin.weight / target_count
        
        
        while True:
            segments = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
            weights = numpy.array(map(operator.attrgetter('weight'), segments))
            cumul_weight = numpy.add.accumulate(weights)
            
            to_merge = segments[cumul_weight <= ideal_weight*self.weight_merge_cutoff]
            if len(to_merge) < 2:
                return
            
            assert cumul_weight[len(to_merge)-1] == sum(weights[:len(to_merge)])
            glom = Segment(weight = cumul_weight[len(to_merge)-1])
            # This selects the proper parent, randomly, according to the relative weights of the
            # particles to be merged. To prove empirically that this works (if logical proof isn't
            # satisfying), do the following:
            #>>> weights = numpy.array([0.01,0.02,0.04,0.2,0.2])
            #>>> weights
            #array([ 0.01,  0.02,  0.04,  0.2 ,  0.2 ])
            #>>> cumul_weight = numpy.add.accumulate(weights)
            #>>> cumul_weight
            #array([ 0.01,  0.03,  0.07,  0.27,  0.47])
            #>>> glom_weight = cumul_weight[3-1]
            #>>> glom_weight
            #0.070000000000000007
            #>>> for i in xrange(0,700000): idx = numpy.digitize((random.uniform(0,glom_weight),),cumul_weight[:3])[0]; counts[idx] += 1;
            #>>> 
            #>>> counts / float(sum(counts))
            #array([ 0.142598,  0.285781,  0.571621])
            #>>> 1/7.0, 2/7.0, 4/7.0
            #(0.14285714285714285, 0.2857142857142857, 0.5714285714285714)
    
            iparent = numpy.digitize((random.uniform(0,glom.weight),),cumul_weight[:len(to_merge)])[0]
            glom.pcoord = to_merge[iparent].pcoord.copy()
            # Set the parent_id field, which tracks trajectory history
            glom.parent_id = to_merge[iparent].parent_id
            # Set the wtg_parent_ids field, which tracks probability flow
            glom.wtg_parent_ids = set()
            for parent in to_merge:
                glom.wtg_parent_ids |= parent.wtg_parent_ids
            
            if log.isEnabledFor(logging.DEBUG):
                log.debug('merging ({:d}) {!r} into {!r}'.format(len(to_merge), to_merge, glom))
            
            for merged in to_merge:
                if merged.initpoint_type == Segment.SEG_INITPOINT_CONTINUES: 
                    self.completed_segs_map[merged.parent_id].endpoint_type = Segment.SEG_ENDPOINT_MERGED
            
            if glom.initpoint_type == Segment.SEG_INITPOINT_CONTINUES:
                self.completed_segs_map[glom.parent_id].endpoint_type = Segment.SEG_ENDPOINT_CONTINUES
                            
            # Remove all particles involved in the merge from the bin
            bin.difference_update(to_merge)
            
            # Add the conglomerate particle resulting from the merge back to the bin
            bin.add(glom)
        
    
    def adjust_count(self, bin):
        '''Adjust the particle count of a bin (by splits/merges) so that it exactly equals the target count'''
        weight_getter = operator.attrgetter('weight')
        
        # split
        while bin.count < bin.target_count:
            # Always split the highest probability particle into two
            particles = sorted(bin, key=weight_getter)
            parent = particles[-1]
            # If the highest-weight particle is newly-created, it will not have a seg_id;
            # In that case, we essentially replace that segment with two new ones
            children = [None] * 2
            for i in xrange(0,2):
                child = Segment(weight = parent.weight/2,
                                pcoord = parent.pcoord.copy(),
                                parent_id = parent.parent_id,
                                wtg_parent_ids = set(parent.wtg_parent_ids))
                children[i] = child              

            if log.isEnabledFor(logging.DEBUG):
                log.debug('splitting {!r} into {:d}: {!r}'.format(parent, len(children), children))             
            bin.remove(parent)
            bin.update(children)
        
        # merge
        while bin.count > bin.target_count:
            # Always merge the two lowest-probability particles
            segments = sorted(bin, key=weight_getter)
            parents = segments[0:2]
            glom_weight = parents[0].weight + parents[1].weight
            parent_ids = [parent.parent_id for parent in parents]
            glom = Segment(weight = glom_weight,
                           wtg_parent_ids = set(parent_ids),
                          )
            
            # Boolean trick: 
            # random(0, glom_weight) < parents[0].weight should map to 0,
            # so just do random(0, glom_weight) >= parents[0].weight
            iparent = (random.uniform(0, glom_weight) >= parents[0].weight)
              
            glom.parent_id = parent_ids[iparent]
            glom.wtg_parent_ids = set()
            for parent in parents:
                glom.wtg_parent_ids |= parent.wtg_parent_ids
            glom.pcoord = parents[iparent].pcoord.copy()
            
            if glom.initpoint_type == Segment.SEG_INITPOINT_CONTINUES:
                self.completed_segs_map[glom.parent_id].endpoint_type = Segment.SEG_ENDPOINT_CONTINUES
                
            # another boolean trick: if index of the selected particle is 0, then the one that gets merged is 1,
            # and vice versa, so just take logical not of iparent as the index
            if parents[~iparent].initpoint_type == Segment.SEG_INITPOINT_CONTINUES:
                self.completed_segs_map[parents[~iparent].parent_id].endpoint_type = Segment.SEG_ENDPOINT_MERGED
            
            if log.isEnabledFor(logging.DEBUG):
                log.debug('merging {:d} {!r} into {!r}'.format(len(parents), parents, glom))
                
            bin.difference_update(parents)
            bin.add(glom)
            
