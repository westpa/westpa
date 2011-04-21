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
    
    def __init__(self, region_set, tfile = None, ltfile = None, edfile = None, fptfile = None,
                 accumulate_statistics = True):
        self.region_set = region_set
        self.tfile = tfile
        self.edfile = edfile
        self.fptfile = fptfile
        self.ltfile = ltfile
        
        self.accumulate_statistics = accumulate_statistics
        self.track_transitions = (tfile is not None)
        self.track_lifetimes = (ltfile is not None)
        self.track_fpts = (fptfile is not None)
        self.track_eds = (edfile is not None)
        
                
        self.n_bins = len(self.region_set.get_all_bins())
        
        self.clear()
        
    def clear(self):
        self.clear_statistics()
        self.clear_state()
        
    def clear_state(self):
        self.last_crossing   = numpy.zeros((self.n_bins,self.n_bins), self.index_dtype)
        self.last_entry      = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_exit       = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_completion = numpy.zeros((self.n_bins,self.n_bins), self.index_dtype)
        self.time_index = None
        self.last_region = None
        
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
                'last_completion': self.last_completion.copy(),
                'time_index':      self.time_index,
                'last_region':     self.last_region}
        
    def set_state(self, state_dict):
        self.last_crossing   = state_dict['last_crossing'].copy()
        self.last_entry      = state_dict['last_entry'].copy()
        self.last_exit       = state_dict['last_exit'].copy()
        self.last_completion = state_dict['last_completion'].copy()
        self.time_index      = state_dict['time_index']
        self.last_region     = state_dict['last_region']
            
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
    
    def accumulate_transitions(self, pcoords, weight = None, time_index=None, continuation = False, label=''):
        """Assign the given progress coordinates to regions, and then determine transitions among regions.
        If continuation is True, ignore the first point (raises an error if accumulate_transitions() has not
        been called at least once already), otherwise use the first point as the initial point."""
        
        if weight is None:
            weights = numpy.ones((len(pcoords),))
        elif numpy.isscalar(weight):
            weights = numpy.zeros((len(pcoords),), dtype=self.weight_dtype)
            weights += weight
        else:
            # weight is an array; check to see if its shape is correct
            if len(weight) != len(pcoords):
                raise ValueError('weight array shape [{!r}] is not compatible with pcoords shape [{!r}]'
                                 .format(weight.shape, pcoords.shape))
            else:
                weights = weight            
        
        indices = self.region_set.map_to_all_indices(pcoords)
        
        # a table of time points, region at that time point, and (yes, believe it), the region at the prior time point
        # in that order
        TIMEPOINT = 0; CURRENT_REGION = 1; LAST_REGION = 2
        
        if continuation:
            tracking_table = numpy.empty((len(pcoords), 3), self.index_dtype)
            tracking_table[0,LAST_REGION] = self.last_region
            tracking_table[:,CURRENT_REGION] = indices
            if len(pcoords) > 1:
                tracking_table[1:,LAST_REGION] = indices[:-1]
            tracking_table[:,TIMEPOINT] = numpy.arange(self.time_index+1, self.time_index+len(pcoords)+1, dtype=self.index_dtype)
        else:
            # not a continuation
            self.clear_state()
            self.time_index = time_index = time_index or 0
            tracking_table = numpy.empty((len(pcoords)-1,3), self.index_dtype)
            tracking_table[:,LAST_REGION] = indices[:-1]
            tracking_table[:,CURRENT_REGION] = indices[1:]
            tracking_table[:,TIMEPOINT] = numpy.arange(self.time_index+1, self.time_index+len(pcoords), dtype=self.index_dtype)
        
        observed_at = tracking_table[:,CURRENT_REGION] != tracking_table[:,LAST_REGION]
        trans_observed = tracking_table[observed_at, :]
        weights_observed = weights[observed_at]

        # pull all references to self out of the inner loops for a 25% reduction in CPU time!
        last_crossing = self.last_crossing
        last_entry = self.last_entry
        last_exit = self.last_exit
        last_completion = self.last_completion
        n_crossings = self.n_crossings
        n_completions = self.n_completions
        record_transition_entry = self.record_transition_entry
        record_lifetime_entry = self.record_lifetime_entry
        record_ed_entry = self.record_ed_entry
        record_fpt_entry = self.record_fpt_entry
        n_bins = self.n_bins
        index_dtype = self.index_dtype
        track_transitions = self.track_transitions
        track_lifetimes = self.track_lifetimes
        track_eds = self.track_eds
        track_fpts = self.track_fpts

        irdisc = numpy.zeros((n_bins,), numpy.bool_)        
        init_regions = numpy.arange(0,n_bins,dtype=index_dtype)
        
        for ((time_index,current_region,last_region), weight) in izip(trans_observed, weights_observed):
            n_crossings[last_region,current_region] += 1
            if track_transitions:
                record_transition_entry(time_index,last_region,current_region,weight, label=label)
            
            if track_lifetimes and last_entry[last_region] > 0:
                lifetime = time_index - last_entry[last_region]
                record_lifetime_entry(time_index, last_region, lifetime, weight, label=label)
                
            irdisc[:] = True
            irdisc[last_region] = False
            irdisc[current_region] = False
            irdisc &= (last_entry > last_completion[:,current_region])
            
            for initial_region in init_regions[irdisc]:
                
                if track_transitions:
                    record_transition_entry(time_index,initial_region,current_region,weight, label=label)
                if track_eds:
                    ed = time_index - last_exit[initial_region]
                    record_ed_entry(time_index,initial_region,current_region,ed,weight, label=label)
                
                if track_fpts and last_completion[current_region,initial_region] > 0:
                        fpt = time_index - last_completion[current_region,initial_region]
                        record_fpt_entry(time_index,initial_region,current_region,fpt,weight, label=label)
                
                last_completion[initial_region,current_region] = time_index
                n_completions[initial_region,current_region] += 1
                         
            last_exit[last_region] = time_index
            last_entry[current_region] = time_index
            last_crossing[last_region,current_region] = time_index
                    
        self.time_index = tracking_table[-1,TIMEPOINT]
        self.last_region = tracking_table[-1,CURRENT_REGION]
