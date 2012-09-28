from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)
import numpy
import operator
from math import ceil
import random
from itertools import izip

import west
from west import Segment
from west.binning import Bin

class ConsistencyError(RuntimeError):
    pass

class AccuracyError(RuntimeError):
    pass

class NewWeightEntry:
    NW_SOURCE_RECYCLED = 0
    
    def __init__(self, source_type, weight, prev_seg_id=None,
                 prev_init_pcoord=None, prev_final_pcoord=None, new_init_pcoord=None,
                 target_state_id = None, initial_state_id = None):
        self.source_type = source_type
        self.weight = weight
        self.prev_seg_id = prev_seg_id
        self.prev_init_pcoord = numpy.asarray(prev_init_pcoord) if prev_init_pcoord is not None else None
        self.prev_final_pcoord = numpy.asarray(prev_final_pcoord) if prev_final_pcoord is not None else None 
        self.new_init_pcoord = numpy.asarray(new_init_pcoord) if new_init_pcoord is not None else None
        self.target_state_id = target_state_id        
        self.initial_state_id = initial_state_id
                

class WEDriver:
    '''A class implemented Huber & Kim's weighted ensemble algorithm over Segment objects.
    This class handles all binning, recycling, and preparation of new Segment objects for the
    next iteration. Binning is accomplished using system.bin_mapper, and per-bin target counts
    are from system.bin_target_counts. 
    
    The workflow is as follows:
    
      1) Call `new_iteration()` every new iteration, providing any recycling targets that are
         in force.
      2) Call `assign()` to assign segments to bins based on their initial and end points. This
         returns the number of walkers that were recycled.
      3) Call `run_we()`, optionally providing a set of initial states that will be used to
         recycle walkers.
         
    Note the presence of flux_matrix, transition_matrix,
    current_iter_segments, next_iter_segments, recycling_segments,
    initial_binning, final_binning, next_iter_binning, and new_weights (to be documented soon).
    '''
    
    weight_split_threshold = 2.0
    weight_merge_cutoff = 1.0
        
    def __init__(self, system=None):
        self.system = system or west.rc.get_system_driver()
                
        self.do_adjust_counts = west.rc.config.get_bool('we.adjust_counts', True)
        log.info('Adjust counts to exactly match target_counts: {}'.format(self.do_adjust_counts))
        
        if 'we.weight_split_threshold' in west.rc.config:
            self.weight_split_threshold = west.rc.config.get_float('we.weight_split_threshold')
        log.info('Split threshold: {}'.format(self.weight_split_threshold))
        
        if 'we.weight_merge_cutoff' in west.rc.config:
            self.weight_merge_cutoff = west.rc.config.get_float('we.weight_merge_cutoff')
        log.info('Merge cutoff: {}'.format(self.weight_merge_cutoff))
        
        # bin mapper and per-bin target counts (see new_iteration for initialization)
        self.bin_mapper = None
        self.bin_target_counts = None
        
        # Target state definitions and corresponding bins
        self.target_states = None
        self.target_state_mask = None
        
        # binning on initial points
        self.initial_binning = None
        
        # binning on final points (pre-WE)
        self.final_binning   = None
                                
        # binning on initial points for next iteration
        self.next_iter_binning = None
                
        # Flux and rate matrices for the current iteration
        self.flux_matrix = None
        self.transition_matrix = None
        
        # Information on new weights (e.g. from recycling) for the next iteration
        self.new_weights = None
        
        # Set of initial states passed to run_we() that are actually used for
        # recycling targets
        self.used_initial_states = None
        
    @property
    def next_iter_segments(self):
        if self.next_iter_binning is None:
            raise RuntimeError('cannot access next iteration segments before running WE')
        
        for _bin in self.next_iter_binning:
            for walker in _bin:
                yield walker
                
    @property
    def current_iter_segments(self):
        for _bin in self.final_binning:
            for walker in _bin:
                yield walker
                
    @property
    def recycling_segments(self):
        if len(self.target_states):
            for target_bin in self.next_iter_binning[self.target_state_mask]:
                for segment in target_bin:
                    yield segment
        else:
            return
                
    def clear(self):
        '''Explicitly delete all Segment-related state.'''
        
        del self.initial_binning, self.final_binning, self.next_iter_binning
        del self.flux_matrix, self.transition_matrix
        del self.new_weights
        
        self.initial_binning = None
        self.final_binning = None
        self.next_iter_binning = None
        self.flux_matrix = None
        self.transition_matrix = None
        self.used_initial_states = None
        self.new_weights = None
        
    def new_iteration(self, target_states=None, new_weights=None):
        '''Prepare for a new iteration. ``initial_states``
        is a sequence of InitialState objects for *available* (i.e. unused) 
        initial states, to be used if any segments result in recycling walkers.
        Target states which generate recycling events are specified in ``target_states``,
        a sequence of TargetState objects. Both ``initial_states`` and ``target_states``
        may be None or empty for equilibrium simulations.
        
        The optional ``new_weight_log`` is a sequence of NewWeightEntry objects which will 
        be used to construct the initial flux matrix.        
        '''
        
        self.clear()
        
        new_weights = new_weights or []
        
        # update mapper, if necessary
        self.bin_mapper = self.system.bin_mapper
        nbins = self.bin_mapper.nbins
        self.bin_target_counts = numpy.array(self.system.bin_target_counts).copy()
        log.debug('mapper is {!r}, handling {:d} bins'.format(self.bin_mapper, nbins))
        
        self.target_states = numpy.array(target_states if target_states else [])
        self.target_state_mask = numpy.zeros((nbins,), numpy.bool_)
        
        self.initial_binning    = self.bin_mapper.construct_bins()
        self.final_binning      = self.bin_mapper.construct_bins()
        self.next_iter_binning  = None
        
        flux_matrix = self.flux_matrix = numpy.zeros((nbins,nbins), dtype=numpy.float64)
        transition_matrix = self.transition_matrix = numpy.zeros((nbins,nbins), numpy.uint)
        
        # map target state specifications to bins
        for tstate in self.target_states:
            tstate_assignment = self.bin_mapper.assign([tstate.pcoord])[0]
            self.target_state_mask[tstate_assignment] = True
            log.debug('target state {!r} mapped to bin {}'.format(tstate, tstate_assignment))

        self.bin_target_counts[self.target_state_mask] = 0            
        
            
        # loop over recycled segments, adding entries to the flux matrix appropriately
        if new_weights:
            init_pcoords = numpy.empty((len(new_weights), self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)
            prev_init_pcoords = numpy.empty((len(new_weights), self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)
            
            for (ientry,entry) in enumerate(new_weights):
                init_pcoords[ientry] = entry.new_init_pcoord
                prev_init_pcoords[ientry] = entry.prev_init_pcoord
            
            init_assignments = self.bin_mapper.assign(init_pcoords)
            prev_init_assignments = self.bin_mapper.assign(prev_init_pcoords)
            
            for (entry, i, j) in izip(new_weights, prev_init_assignments, init_assignments):
                flux_matrix[i,j] += entry.weight
                transition_matrix[i,j] += 1
                
        self.recycling_log = None
                        
    def assign(self, segments, initializing=False):
        '''Assign segments to initial and final bins.
        Returns the number of initial states that must be generated.
        If ``initializing`` is True, then the "final" bin assignments will be identical to the initial bin
        assignments, a condition required for seeding a new iteration from pre-existing segments.'''

        # collect initial and final coordinates into one place        
        all_pcoords = numpy.empty((2,len(segments), self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)
        
        for iseg, segment in enumerate(segments):
            all_pcoords[0,iseg] = segment.pcoord[0,:]
            all_pcoords[1,iseg] = segment.pcoord[-1,:]
        
        # assign based on initial and final progress coordinates
        initial_assignments = self.bin_mapper.assign(all_pcoords[0,:,:])
        if initializing:
            final_assignments = initial_assignments
        else:
            final_assignments = self.bin_mapper.assign(all_pcoords[1,:,:])
            
        initial_binning = self.initial_binning
        final_binning = self.final_binning
        flux_matrix = self.flux_matrix
        transition_matrix = self.transition_matrix
        for (segment,iidx,fidx) in izip(segments, initial_assignments, final_assignments):
            initial_binning[iidx].add(segment)
            final_binning[fidx].add(segment)
            flux_matrix[iidx,fidx] += segment.weight
            transition_matrix[iidx,fidx] += 1
            
        n_recycled_total = sum(len(_bin) for _bin in self.final_binning[self.target_state_mask])
        return n_recycled_total
    
    def _recycle_walkers(self, initial_states):
        '''Recycle walkers'''
        
        # recall that every walker we deal with is already a new segment, so to recycle, we actually move
        # the appropriate Segment from the target bin to the initial state bin
        
        self.new_weights = []
        self._avail_initial_states = {state.state_id: state for state in initial_states or []}
        self._used_initial_states = {}
        
        n_recycled_walkers = len(list(self.recycling_segments))
        if not n_recycled_walkers:
            return
        elif n_recycled_walkers > len(initial_states):
            raise ConsistencyError('need {} initial states for recycling, but only {} present'
                                   .format(n_recycled_walkers,len(initial_states)))

        istateiter = iter(initial_states)
        for (itarget, target_bin) in enumerate(self.next_iter_binning[self.target_state_mask]):
            for segment in set(target_bin):
                target_state = self.target_states[itarget]
                initial_state = istateiter.next()
                if log.isEnabledFor(logging.DEBUG):
                    log.debug('recycling {!r} to initial state {!r}'.format(segment, initial_state))
                
                istate_assignment = self.bin_mapper.assign([initial_state.pcoord])[0]
                parent = self.parent_map[segment.parent_id]
                parent.endpoint_type = Segment.SEG_ENDPOINT_RECYCLED
                segment.parent_id = -(initial_state.state_id+1)
                segment.prev_init_pcoord = parent.pcoord[0]
                segment.pcoord[0] = initial_state.pcoord

                self.new_weights.append(NewWeightEntry(source_type=NewWeightEntry.NW_SOURCE_RECYCLED,
                                                       weight=parent.weight, prev_seg_id=parent.seg_id, 
                                                       prev_init_pcoord=parent.pcoord[0],
                                                       prev_final_pcoord=parent.pcoord[-1],
                                                       new_init_pcoord=initial_state.pcoord,
                                                       target_state_id=target_state.state_id,
                                                       initial_state_id=initial_state.state_id) )


                self.next_iter_binning[istate_assignment].add(segment)
                self._used_initial_states[initial_state.state_id] = initial_state
                initial_state.iter_used = segment.n_iter
                log.debug('marking initial state {!r} as used'.format(initial_state))
                del self._avail_initial_states[initial_state.state_id]
                target_bin.remove(segment)                
                 
            assert len(target_bin) == 0
                            
    def _split_walker(self, segment, m, bin):
        '''Split the walker ``segment`` (in ``bin``) into ``m`` walkers'''

        bin.remove(segment)
            
        new_segments = []
        for _inew in xrange(0,m):
            new_segment = Segment(n_iter = segment.n_iter, #previously incremented
                                  weight = segment.weight/m,
                                  parent_id = segment.parent_id,
                                  wtg_parent_ids = set(segment.wtg_parent_ids),
                                  pcoord = segment.pcoord.copy(),
                                  status = Segment.SEG_STATUS_PREPARED)
            new_segment.pcoord[0,:] = segment.pcoord[0,:]
            new_segments.append(new_segment)
            
        bin.update(new_segments)

        if log.isEnabledFor(logging.DEBUG):
            log.debug('splitting {!r} into {:d}:\n    {!r}'.format(segment, m, new_segments))
            
        return new_segments 
                        
    def _merge_walkers(self, segments, cumul_weight, bin):
        '''Merge the given ``segments`` in ``bin``, previously sorted by weight, into one conglomerate segment.
        ``cumul_weight`` is the cumulative sum of the weights of the ``segments``; this may be None to calculate here.'''
        
        if cumul_weight is None:
            cumul_weight = numpy.add.accumulate([segment.weight for segment in segments])
        
        glom = Segment(n_iter = segments[0].n_iter, # assumed correct (and equal among all segments)
                       weight = cumul_weight[len(segments)-1],
                       status = Segment.SEG_STATUS_PREPARED,
                       pcoord = self.system.new_pcoord_array(),
                       )
                
        # Select the history to use
        # The following takes a random number in the interval 0 <= x < glom.weight, then
        # sees where this value falls among the (sorted) weights of the segments being merged;
        # this ensures that a walker with (e.g.) twice the weight of its brethren has twice the
        # probability of having its history selected for continuation
        iparent = numpy.digitize((random.uniform(0,glom.weight),),cumul_weight)[0]
        gparent_seg = segments[iparent]
        
        # Inherit history from this segment ("gparent" stands for "glom parent", as opposed to historical
        # parent). 
        glom.parent_id = gparent_seg.parent_id
        glom.pcoord[0,:] = gparent_seg.pcoord[0,:]
        
        # Weight comes from all segments being merged, and therefore all their
        # parent segments
        glom.wtg_parent_ids = set()
        for segment in segments:
            glom.wtg_parent_ids |= segment.wtg_parent_ids
            
        # The historical parent of gparent is continued; all others are marked as merged
        for segment in segments:
            if segment is gparent_seg:
                # we must ignore initial states here...
                if segment.parent_id >= 0:
                    self.parent_map[segment.parent_id].endpoint_type = Segment.SEG_ENDPOINT_CONTINUES
            else:
                # and "unuse" an initial state here
                if segment.parent_id >= 0:
                    self.parent_map[segment.parent_id].endpoint_type = Segment.SEG_ENDPOINT_MERGED
                else:
                    initial_state = self._used_initial_states.pop(segment.initial_state_id)
                    log.debug('freeing initial state {!r} for future use (merged)'.format(initial_state))
                    self._avail_initial_states[initial_state.state_id] = initial_state
                    initial_state.iter_used = None

        if log.isEnabledFor(logging.DEBUG):
            log.debug('merging ({:d}) {!r} into 1:\n    {!r}'.format(len(segments), segments, glom))
                
        bin.difference_update(segments)
        bin.add(glom)
        
    def _split_by_weight(self, ibin):
        '''Split overweight particles'''
        
        bin = self.next_iter_binning[ibin]
        target_count = self.bin_target_counts[ibin]
        segments = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
        weights = numpy.array(map(operator.attrgetter('weight'), segments))
        ideal_weight = weights.sum() / target_count
 
        if len(bin) > 0:
            assert target_count > 0
               
        to_split = segments[weights > self.weight_split_threshold*ideal_weight]
        
        for segment in to_split:
            m = int(ceil(segment.weight / ideal_weight))
            self._split_walker(segment, m, bin)

    
    def _merge_by_weight(self, ibin):
        '''Merge underweight particles'''
        
        bin = self.next_iter_binning[ibin]
        target_count = self.bin_target_counts[ibin]
        weight = sum(map(operator.attrgetter('weight'), bin))
        target_count = self.bin_target_counts[ibin]
        ideal_weight = weight / target_count    
        
        while True:
            segments = numpy.array(sorted(bin, key=operator.attrgetter('weight')), dtype=numpy.object_)
            weights = numpy.array(map(operator.attrgetter('weight'), segments))
            cumul_weight = numpy.add.accumulate(weights)
            
            to_merge = segments[cumul_weight <= ideal_weight*self.weight_merge_cutoff]
            if len(to_merge) < 2:
                return
            
            self._merge_walkers(to_merge, cumul_weight, bin)
    
    def _adjust_count(self, ibin):
        bin = self.next_iter_binning[ibin]
        target_count = self.bin_target_counts[ibin]
        weight_getter = operator.attrgetter('weight')

        # split        
        while len(bin) < target_count:
            log.debug('adjusting counts by splitting')
            # always split the highest probability walker into two
            segments = sorted(bin, key=weight_getter)
            self._split_walker(segments[-1], 2, bin)
            
        # merge
        while len(bin) > target_count:
            log.debug('adjusting counts by merging')
            # always merge the two lowest-probability walkers
            segments = sorted(bin, key=weight_getter)
            self._merge_walkers(segments[:2], cumul_weight=None, bin=bin)

    def _check_pre(self):
        for ibin, _bin in enumerate(self.next_iter_binning):
            if self.bin_target_counts[ibin] == 0 and len(_bin) > 0:
                raise ConsistencyError('bin {:d} has target count of 0 but contains {:d} walkers'.format(ibin, len(_bin)))

    def _check_post(self):
        for segment in self.next_iter_segments:
            if segment.weight == 0:
                raise ConsistencyError('segment {!r} has weight of zero')
                                    
    def run_we(self, initial_states=None):
        '''Run weighted ensemble recycling and split/merge on the segments previously assigned to
        bins using ``assign_segments``. Enough unused initial states must be present in 
        ``self.initial_states`` for every recycled walker to be assigned an initial state.
        After this function completes, ``self.flux_matrix`` contains a valid flux matrix for this
        iteration (including any contributions from recycling from the previous iteration), and
        ``self.next_iter_segments`` contains a list of segments ready for the next iteration,
        with appropriate values set for weight, endpoint type, parent walkers, and so on.'''
        
        self.parent_map = {}
        self.next_iter_binning = self.bin_mapper.construct_bins()

        # Create new segments for the next iteration        
        # We assume that everything is going to continue without being touched by recycling or WE, and
        # adjust later
        new_pcoord_array = self.system.new_pcoord_array
        n_iter = None
        for ibin, _bin in enumerate(self.final_binning):
            for segment in _bin:
                if n_iter is None:
                    n_iter = segment.n_iter
                else:
                    assert segment.n_iter == n_iter
                    
                segment.endpoint_type = Segment.SEG_ENDPOINT_CONTINUES
                new_segment = Segment(n_iter=segment.n_iter+1,
                                      parent_id=segment.seg_id,
                                      weight=segment.weight,
                                      wtg_parent_ids=[segment.seg_id],
                                      pcoord=new_pcoord_array(),
                                      status=Segment.SEG_STATUS_PREPARED)
                new_segment.pcoord[0] = segment.pcoord[-1]
                self.next_iter_binning[ibin].add(new_segment)
                
                # Store a link to the parent segment, so we can update its endpoint status as we need,
                # based on its ID
                self.parent_map[segment.seg_id] = segment
                                
        
        # Then, recycle walkers that exist in target states 
        self._recycle_walkers(initial_states)
        
        # sanity check
        self._check_pre()
        
        # Regardless of current particle count, always split overweight particles and merge underweight particles
        # Then and only then adjust for correct particle count
        for (ibin,bin) in enumerate(self.next_iter_binning):
            if len(bin) == 0: 
                continue
                        
            self._split_by_weight(ibin)
            self._merge_by_weight(ibin)
            if self.do_adjust_counts:
                self._adjust_count(ibin)
            
        self._check_post()
        
        self.new_weights = self.new_weights or []
        self.used_initial_states = set(self._used_initial_states.itervalues())
        self.avail_initial_states = set(self._avail_initial_states.itervalues())
        
        log.debug('used initial states: {!r}'.format(self.used_initial_states))
        log.debug('available initial states: {!r}'.format(self.avail_initial_states))

        
    
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
                    
            
