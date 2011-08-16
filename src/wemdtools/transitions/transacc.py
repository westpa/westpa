from __future__ import division; __metaclass__ = type
import numpy
from itertools import izip

import logging
log = logging.getLogger(__name__)

from wemdtools.stats.accumulator import RunningStatsAccumulator

# Indices for fields in transition data tuples
TDAT_NITER = 0          # WE iteration in which the transition occurred
TDAT_TIMEPOINT = 1      # Time within the trajectory 
TDAT_INIT_REGION = 2    # Index of initial region of transition 
TDAT_FINAL_REGION = 3   # Index of final region of transition
TDAT_WEIGHT = 4         # Weight at time of crossing into final region
TDAT_IPROB = 5          # Weight at time of crossing out of initial region
TDAT_LIFETIME = 6       # Lifetime (of being in TDAT_INIT_REGION) ended by this crossing
TDAT_ED = 7             # Event duration of this transition              

# Numpy data type
tdat_dtype = numpy.dtype( [('n_iter', numpy.int32),
                           ('timepoint', numpy.uint64),
                           ('init_region', numpy.uint16),
                           ('final_region', numpy.uint16),
                           ('weight', numpy.float64),
                           ('iprob', numpy.float64),
                           ('lifetime', numpy.uint64),
                           ('ed', numpy.uint64),
                          ])

class TransitionEventAccumulator:
    index_dtype  = numpy.uintp
    count_dtype  = numpy.uint64
    time_dtype   = numpy.float64
    weight_dtype = numpy.float64
    output_tdat_chunksize = 10000    # HDF5 chunksize for transition data (~500 KiB)
    tdat_buffersize = 1000000         # Internal buffer length (~50 MiB)
    
    def __init__(self, region_set, output_group):
        self.region_set = region_set
        self.n_bins = len(self.region_set.get_all_bins())
        self.init_regions = numpy.arange(0,self.n_bins, dtype=self.index_dtype)
        
        # HDF5 group in which to store results
        self.output_group = output_group
        self.tdat_buffer = numpy.empty((self.tdat_buffersize,), dtype=tdat_dtype)
        self.tdat_buffer_offset = 0
        self.output_tdat_offset = 0
        self.output_tdat_ds = None
                        
        # Accumulators/counters
        self.n_crossings     = None # shape (n_bins,n_bins)
        self.n_completions   = None # (n_bins,n_bins)        
        self.lt_acc          = None # (n_bins,)
        self.ed_acc          = None # (n_bins, n_bins)
        self.fpt_acc         = None # (n_bins,n_bins)

        # Time points and per-timepoint data
        self.last_crossing   = None # shape (n_bins,n_bins)
        self.last_entry      = None # (n_bins,)
        self.last_exit       = None # (n_bins,)
        self.last_exit_weights = None # (n_bins,)        
        self.last_completion = None # (n_bins,n_bins)
        
        # Analysis continuation information
        self.time_index         = None # current time index for separate calls on same trajectory
        self.last_region        = None # last region occupied, for separate calls on same trajectory
        self.last_region_weight = None # total weight in self.last_region at end of last processing step
        
        self.clear()
        
    def clear(self):
        self.clear_statistics()
        self.clear_state()
        self.tdat_buffer = numpy.empty((self.tdat_buffersize,), dtype=tdat_dtype)
        self.tdat_buffer_offset = 0
        self.output_tdat_offset = 0
        self.output_tdat_ds = None
        
    def clear_state(self):
        self.last_crossing      = numpy.zeros((self.n_bins,self.n_bins), self.index_dtype)
        self.last_entry         = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_exit          = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_exit_weights  = numpy.zeros((self.n_bins,), self.weight_dtype)
        self.last_completion    = numpy.zeros((self.n_bins,self.n_bins), self.index_dtype)
        self.time_index         = None
        self.last_region        = None
        self.last_region_weight = None
        
    def clear_statistics(self):
        self.n_crossings     = numpy.zeros((self.n_bins,self.n_bins), self.count_dtype)
        self.n_completions   = numpy.zeros((self.n_bins,self.n_bins), self.count_dtype)        
        self.lt_acc  = RunningStatsAccumulator(shape=(self.n_bins,), dtype=self.time_dtype, count_dtype=self.count_dtype,
                                               weight_dtype = self.weight_dtype, mask_value=0)
        self.ed_acc  = RunningStatsAccumulator(shape=(self.n_bins,self.n_bins), dtype=self.time_dtype, count_dtype=self.count_dtype,
                                               weight_dtype = self.weight_dtype, mask_value=0)
        self.fpt_acc = RunningStatsAccumulator(shape=(self.n_bins,self.n_bins), dtype=self.time_dtype, count_dtype=self.count_dtype,
                                               weight_dtype = self.weight_dtype, mask_value=0)

    def get_state(self):
        return {'last_crossing':   self.last_crossing.copy(),
                'last_entry':      self.last_entry.copy(),
                'last_exit':       self.last_exit.copy(),
                'last_exit_weights': self.last_exit_weights.copy(),
                'last_completion': self.last_completion.copy(),
                'time_index':      self.time_index,
                'last_region':     self.last_region,
                'last_region_weight': self.last_region_weight}
        
    def set_state(self, state_dict):
        self.last_crossing   = state_dict['last_crossing'].copy()
        self.last_entry      = state_dict['last_entry'].copy()
        self.last_exit       = state_dict['last_exit'].copy()
        self.last_exit_weights = state_dict['last_exit_weights'].copy()
        self.last_completion = state_dict['last_completion'].copy()
        self.time_index      = state_dict['time_index']
        self.last_region     = state_dict['last_region']
        self.last_region_weight = state_dict['last_region_weight']
        
        
    def record_transition_data(self, tdat):
        """Update running statistics and write transition data to HDF5 (with buffering)"""
        
        # Update running statistics
        # event durations are weighted, everything else is unweighted
        lt_count = self.lt_acc.count
        lt_weight = self.lt_acc.weight
        lt_sum = self.lt_acc.sum
        lt_sqsum = self.lt_acc.sqsum
        
        ed_count = self.ed_acc.count
        ed_weight = self.ed_acc.weight
        ed_sum = self.ed_acc.sum
        ed_sqsum = self.ed_acc.sqsum
                
        for (n_iter, timepoint, init_region, final_region, weight, iprob, lifetime, ed) in tdat:
            idx = (init_region, final_region)
                                    
            if lifetime > 0:
                lt_count[init_region] += 1
                lt_weight[init_region] += 1
                lt_sum[init_region] += weight*lifetime
                lt_sqsum[init_region] += weight*lifetime*lifetime
            
            if ed > 1:
                ed_count[idx] += 1
                ed_weight[idx] += weight
                ed_sum[idx] += weight*ed
                ed_sqsum[idx] += weight*ed*ed
        
        
        # Write out accumulated transition data
        if self.output_tdat_ds is None:        
            # Create dataset
            self.output_tdat_ds = self.output_group.create_dataset('transitions', shape=(self.output_tdat_chunksize,), 
                                                                   dtype=tdat_dtype, maxshape=(None,), 
                                                                   chunks=(self.output_tdat_chunksize,))
            
        # If the amount of data to write exceeds our remaining buffer space, flush the buffer, then
        # write data directly to HDF5, otherwise just add to the buffer and wait for the last flush
        if len(tdat) + self.tdat_buffer_offset > self.tdat_buffersize:
            self.flush_transition_data()
            ub = self.output_tdat_offset + len(tdat)
            self.output_tdat_ds.resize((ub,))
            self.output_tdat_ds[self.output_tdat_offset:ub] = tdat
            self.output_tdat_offset += len(tdat)
        else:
            self.tdat_buffer[self.tdat_buffer_offset:(self.tdat_buffer_offset+len(tdat))] = tdat
            self.tdat_buffer_offset += len(tdat)
        
    def flush_transition_data(self):
        """Flush any unwritten output that may be present"""
        if self.output_tdat_ds is None:
            return
                
        # self.tdat_buffer_offset is the number of items in the buffer
        nbuf = self.tdat_buffer_offset
        if nbuf == 0: return
        ub = nbuf + self.output_tdat_offset
        if ub > self.output_tdat_ds.len():
            # Resize dataset to fit data
            self.output_tdat_ds.resize((ub,))
            
        self.output_tdat_ds[self.output_tdat_offset:ub] = self.tdat_buffer[:nbuf]
        self.output_tdat_offset += nbuf
        self.tdat_buffer_offset = 0
        
                    
    def accumulate_transitions(self, pcoords, 
                                   weight = None, region_weights = None, 
                                   time_index=None, continuation = False,
                                   n_iter=None):
        """Assign the given progress coordinates to regions, and then determine transitions among regions.
        If continuation is True, ignore the first point (raises an error if accumulate_transitions() has not
        been called at least once already), otherwise use the first point as the initial point."""
        
        # Coerce the given weight into a vector the same length as pcoords
        if weight is None:
            weights = numpy.ones((len(pcoords),), dtype=self.weight_dtype)
        elif numpy.isscalar(weight):
            weights = numpy.empty((len(pcoords),), dtype=self.weight_dtype)
            weights[...] = weight
        else:
            # weight is an array; check to see if its shape is correct
            if len(weight) != len(pcoords):
                raise ValueError('weight array shape [{!r}] is not compatible with pcoords shape [{!r}]'
                                 .format(weight.shape, pcoords.shape))
            else:
                weights = weight
        
        if region_weights is None:
            region_weights = numpy.ones((len(pcoords), self.n_bins))
                
        indices = self.region_set.map_to_all_indices(pcoords)
        
        # a table of time points, region at that time point, and (yes, believe it), the region at the prior time point
        # in that order; the idea is to do as much vectorized work as possible, even down to the comparison of whether
        # we have progressed from one region to another from one time point to the next
        tracking_table = None
        TIMEPOINT = 0; CURRENT_REGION = 1; LAST_REGION = 2
        
        # The total weight in LAST_REGION at the previous time point, used for calculating rates
        region_weight_table = None
        
        if continuation:
            # continuation; take data from end of prior call
            tracking_table = numpy.empty((len(pcoords), 3), self.index_dtype)
            tracking_table[0,LAST_REGION] = self.last_region
            tracking_table[:,CURRENT_REGION] = indices
            region_weight_table = numpy.empty((len(pcoords),), self.weight_dtype)
            region_weight_table[0] = self.last_region_weight
            if len(pcoords) > 1:
                tracking_table[1:,LAST_REGION] = indices[:-1]
                #region_weight_table[1:] = region_weights[:-1,indices[:-1]]
                region_weight_table[1:] = [region_weights[i,indices[i]] for i in xrange(0,len(indices)-1)]
            til = long(self.time_index)
            xr = xrange(til+1, til+len(pcoords)+1)
            tracking_table[:,TIMEPOINT] = xr
        else:
            # not a continuation; initial point comes from given data, not from end point of prior call
            self.clear_state()
            self.time_index = time_index = time_index or 0
            tracking_table = numpy.empty((len(pcoords)-1,3), self.index_dtype)
            tracking_table[:,LAST_REGION] = indices[:-1]
            tracking_table[:,CURRENT_REGION] = indices[1:]
            region_weight_table = numpy.empty((len(pcoords)-1,), self.weight_dtype)
            region_weight_table[:] = [region_weights[i,indices[i]] for i in xrange(0,len(indices)-1)]                        
            til = long(self.time_index)
            xr = xrange(til+1, til+len(pcoords))
            tracking_table[:,TIMEPOINT] = xr
        
        observed_at = tracking_table[:,CURRENT_REGION] != tracking_table[:,LAST_REGION]
        trans_observed = tracking_table[observed_at, :]
        weights_observed = weights[observed_at]
        last_region_weights_observed = region_weight_table[observed_at]

        # pull all references to self out of the inner loops for a 25% reduction in CPU time (!)
        last_crossing = self.last_crossing
        last_entry = self.last_entry
        last_exit = self.last_exit
        last_exit_weights = self.last_exit_weights
        last_completion = self.last_completion
        n_crossings = self.n_crossings
        n_completions = self.n_completions
        n_bins = self.n_bins
        n_iter = n_iter or 0
                
        tdat = []

        # a discriminator indicating which initial regions we have to deal with when calculating transitions
        irdisc = numpy.zeros((n_bins,), numpy.bool_)        
        init_regions = self.init_regions
        
        for ((time_index,current_region,last_region), weight, last_region_weight) in izip(trans_observed, weights_observed, 
                                                                                          last_region_weights_observed):
            n_crossings[last_region,current_region] += 1
            n_completions[last_region,current_region] += 1

            if last_entry[last_region] > 0:
                lifetime = time_index - last_entry[last_region]
            else:
                lifetime = 0
                
            tdat.append((n_iter, time_index, last_region, current_region, weight, last_region_weight, lifetime, 1))
                        
            # Note that transitions, event durations, etc. all require a transition *through* an intermediate state
            # i.e. the transitions file and count matrix do not reflect transitions across one (or more) boundaries
            # with no intervening time (thanks to Brandon Mills for pointing out how confusing this is without
            # a few more comments)                
            irdisc[:] = True # We need to look at everything
            #irdisc[last_region] = False # except transitions beginning in the last region (that's just a crossing, handled above)
            irdisc &= (last_entry > last_completion[:,current_region]) # and only those visited since the last time we were here
            
            for initial_region in init_regions[irdisc]:
                if last_exit[initial_region] > 0:        
                    ed = time_index - last_exit[initial_region]                        
                    tdat.append((n_iter, time_index, initial_region, current_region, weight, last_exit_weights[initial_region], 0, ed))
                
                last_completion[initial_region,current_region] = time_index
                if initial_region != last_region:
                    n_completions[initial_region,current_region] += 1
                         
            last_exit[last_region] = time_index
            last_exit_weights[last_region] = last_region_weight
            last_entry[current_region] = time_index
            last_crossing[last_region,current_region] = time_index
                    
        # Update the last-point record for later continuation
        self.time_index = tracking_table[-1,TIMEPOINT]
        self.last_region = tracking_table[-1,CURRENT_REGION]
        self.last_region_weight = region_weights[-1, self.last_region]
        
        if tdat:
            self.record_transition_data(tdat)
            

