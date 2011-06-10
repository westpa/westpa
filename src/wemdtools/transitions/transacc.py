from __future__ import division; __metaclass__ = type
import numpy
from itertools import izip

import logging
log = logging.getLogger(__name__)

from wemdtools.stats.accumulator import RunningStatsAccumulator

class TransitionEventAccumulator:
    index_dtype  = numpy.uintp
    count_dtype  = numpy.uint64
    time_dtype   = numpy.float64
    weight_dtype = numpy.float64
    
    def __init__(self, region_set, tfile = None, ltfile = None, edfile = None, fptfile = None, ratefile = None,
                 accumulate_statistics = True):
        self.region_set = region_set
        self.n_bins = len(self.region_set.get_all_bins())
        self.init_regions = numpy.arange(0,self.n_bins, dtype=self.index_dtype)
        
        self.tfile = tfile       # file to which to write an accounting of transitions
        self.edfile = edfile     # event duration file
        self.fptfile = fptfile   # first passage time file
        self.ltfile = ltfile     # lifetime (dwell time) file
        self.ratefile = ratefile # bin-to-bin rates 
        
        # whether or not to record various kinds of data -- each takes a little bit of performance out of 
        # the overall procedure
        self.accumulate_statistics = accumulate_statistics # whether or not to accumulate running averages/std.devs.
        self.track_transitions = (tfile is not None)       # whether to record individual transitions
        self.track_lifetimes = (ltfile is not None)        # whether to record lifetimes
        self.track_fpts = (fptfile is not None)            # whether to record first passage times
        self.track_eds = (edfile is not None)              # whether to record event duration times
        self.track_rates = (ratefile is not None and edfile is not None) # whether to record rates (requires tracking eds)
        
        # Accumulators/counters
        self.n_crossings     = None # shape (n_bins,n_bins)
        self.n_completions   = None # (n_bins,n_bins)        
        self.lt_acc          = None # (n_bins,)
        self.ed_acc          = None # (n_bins, n_bins)
        self.fpt_acc         = None # (n_bins,n_bins)
        self.rate_acc        = None # (n_bins,n_bins)
        self.flux_acc        = None # (n_bins, n_bins)

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
        self.rate_acc = RunningStatsAccumulator(shape=(self.n_bins,self.n_bins), 
                                                dtype=self.weight_dtype, count_dtype=self.count_dtype,
                                                weight_dtype = self.weight_dtype, mask_value = 0)
        self.flux_acc = RunningStatsAccumulator(shape=(self.n_bins,self.n_bins),
                                                dtype=self.weight_dtype, count_dtype=self.count_dtype,
                                                weight_dtype = self.weight_dtype, mask_value = 0)

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
            
    def record_transition_entry(self, time_index, initial_region, final_region, weight, label=''):
        if self.tfile:
            tlabel = '{:d}->{:d}'.format(int(initial_region),int(final_region))
            self.tfile.write('{label:<24s}    {time_index:20d}    {tlabel:<24s}    {weight:20.14e}\n'
                             .format(time_index=long(time_index), tlabel=tlabel, weight=float(weight), label=label))
    
    def record_lifetime_entry(self, time_index, final_region, lifetime, weight, label=''):
        if self.accumulate_statistics:
            lt_acc = self.lt_acc
            lt_acc.count[final_region] += 1
            lt_acc.weight[final_region] += weight
            lt_acc.sum[final_region] += weight*lifetime
            lt_acc.sqsum[final_region] += weight*lifetime*lifetime
        
        if self.ltfile:
            self.ltfile.write('{label:<24s}    {time_index:20d}    {slabel:10d}    {lifetime:20d}    {weight:20.14e}\n'
                              .format(time_index=long(time_index), slabel=long(final_region), lifetime=long(lifetime), 
                                      weight=float(weight), label=label))
    
    def record_ed_entry(self, time_index, initial_region, final_region, ed, weight, label=''):
        if self.accumulate_statistics:
            ed_acc = self.ed_acc
            ed_acc.count[initial_region,final_region] += 1
            ed_acc.weight[initial_region,final_region] += weight
            ed_acc.sum[initial_region,final_region] += weight*ed
            ed_acc.sqsum[initial_region,final_region] += weight*ed*ed
        if self.edfile:
            tlabel = '{:d}->{:d}'.format(int(initial_region),int(final_region))
            self.edfile.write('{label:<24s}    {time_index:20d}    {tlabel:<24s}    {ed:20d}    {weight:20.14e}\n'
                              .format(time_index=long(time_index), tlabel=tlabel, ed=long(ed), weight=float(weight),
                                      label=label))
    
    def record_fpt_entry(self, time_index, initial_region, final_region, fpt, weight, label=''):
        if self.accumulate_statistics:
            fpt_acc = self.fpt_acc
            fpt_acc.count[initial_region,final_region] += 1
            fpt_acc.weight[initial_region,final_region] += weight
            fpt_acc.sum[initial_region,final_region] += weight*fpt
            fpt_acc.sqsum[initial_region,final_region] += weight*fpt*fpt
        if self.fptfile:
            tlabel = '{:d}->{:d}'.format(int(initial_region),int(final_region))
            self.fptfile.write('{label:<24s}    {time_index:20d}    {tlabel:<24s}    {fpt:20d}    {weight:20.14e}\n'
                              .format(time_index=long(time_index), tlabel=tlabel, fpt=long(fpt), weight=float(weight), 
                                      label=label))
            
    def record_rate_entry(self, time_index, initial_region, final_region, flux, dt, pop, rate, label=''):
        if self.accumulate_statistics:
            rate_acc = self.rate_acc
            rate_acc.count[initial_region,final_region] +=1
            rate_acc.weight[initial_region,final_region] += 1.0
            rate_acc.sum[initial_region,final_region] += rate
            rate_acc.sqsum[initial_region,final_region] += rate*rate
            
            flux_acc = self.flux_acc
            flux_acc.count[initial_region,final_region] += 1
            flux_acc.weight[initial_region,final_region] += 1.0
            flux_acc.sum[initial_region,final_region] += flux
            flux_acc.sqsum[initial_region,final_region] += flux*flux
            
        if self.ratefile:
            tlabel = '{:d}->{:d}'.format(int(initial_region), int(final_region))
            self.ratefile.write('{label:<24s}    {time_index:20d}    {tlabel:<24s}    {flux:20.14e}    {dt:20d}    {pop:20.14e}    {rate:20.14e}\n'
                                .format(time_index=long(time_index), tlabel=tlabel, 
                                        flux=float(flux), dt=long(dt), pop=float(pop), rate=float(rate), label=label))

    def accumulate_transitions(self, pcoords, 
                                   weight = None, region_weights = None, 
                                   time_index=None, continuation = False, label=''):
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
        record_transition_entry = self.record_transition_entry
        record_lifetime_entry = self.record_lifetime_entry
        record_ed_entry = self.record_ed_entry
        record_fpt_entry = self.record_fpt_entry
        record_rate_entry = self.record_rate_entry
        n_bins = self.n_bins
        track_transitions = self.track_transitions
        track_lifetimes = self.track_lifetimes
        track_eds = self.track_eds
        track_fpts = self.track_fpts
        track_rates = self.track_rates

        # a discriminator indicating which initial regions we have to deal with when calculating transitions
        irdisc = numpy.zeros((n_bins,), numpy.bool_)        
        init_regions = self.init_regions
        
        for ((time_index,current_region,last_region), weight, last_region_weight) in izip(trans_observed, weights_observed, 
                                                                                          last_region_weights_observed):
            n_crossings[last_region,current_region] += 1
            if track_transitions:
                record_transition_entry(time_index,last_region,current_region,weight, label=label)
            
            if track_lifetimes and last_entry[last_region] > 0:
                lifetime = time_index - last_entry[last_region]
                record_lifetime_entry(time_index, last_region, lifetime, weight, label=label)
                
            if track_rates and last_region_weight > 0:
                # Fluxes for crossings are always divided by a time of one (the resolution of pcoord data)
                rate = weight / last_region_weight
                #record_rate_entry(time_index, last_region, current_region, rate, label=label)
                #record_rate_entry(self, time_index, initial_region, final_region, flux, dt, pop, rate, label='')
                record_rate_entry(time_index, last_region, current_region, flux=weight, dt=1, pop=last_region_weight,
                                  rate=rate, label=label)
            
            # Note that transitions, event durations, etc. all require a transition *through* an intermediate state
            # i.e. the transitions file and count matrix do not reflect transitions across one (or more) boundaries
            # with no intervening time (thanks to Brandon Mills for pointing out how confusing this is without
            # a few more comments)                
            irdisc[:] = True # We need to look at everything
            irdisc[last_region] = False # except transitions beginning in the last region (that's just a crossing, reported above)
            #irdisc[current_region] = False # and transitions beginning in the current region (that's coming full circle)
            irdisc &= (last_entry > last_completion[:,current_region]) # and only those visited since the last time we were here
            
            for initial_region in init_regions[irdisc]:
                if track_transitions:
                    record_transition_entry(time_index,initial_region,current_region,weight, label=label)
                if track_eds:
                    ed = time_index - last_exit[initial_region]
                    record_ed_entry(time_index,initial_region,current_region,ed,weight, label=label)
                if track_fpts and last_completion[current_region,initial_region] > 0:
                    fpt = time_index - last_completion[current_region,initial_region]
                    record_fpt_entry(time_index,initial_region,current_region,fpt,weight, label=label)
                if track_rates and track_eds and last_exit[initial_region] > 0 and last_exit_weights[initial_region] > 0:
                    # k_ij = 1/(t_j - t_i) * (Prob arriving in j) / (Prob of leaving i)
                    rate = weight / last_exit_weights[initial_region] / ed
                    #record_rate_entry(time_index, initial_region, current_region, rate, label=label)
                    record_rate_entry(time_index, initial_region, current_region, flux=weight, dt=ed, 
                                      pop=last_exit_weights[initial_region], rate=rate, label=label)
                
                last_completion[initial_region,current_region] = time_index
                n_completions[initial_region,current_region] += 1
                         
            last_exit[last_region] = time_index
            last_exit_weights[last_region] = last_region_weight
            last_entry[current_region] = time_index
            last_crossing[last_region,current_region] = time_index
                    
        # Update the last-point record for later continuation
        self.time_index = tracking_table[-1,TIMEPOINT]
        self.last_region = tracking_table[-1,CURRENT_REGION]
        self.last_region_weight = region_weights[-1, self.last_region]
