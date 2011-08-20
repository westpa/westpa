from __future__ import division; __metaclass__ = type
import numpy
from itertools import izip

import logging
log = logging.getLogger(__name__)

class TransitionEventAccumulator:
    index_dtype  = numpy.uintp
    count_dtype  = numpy.uint64
    weight_dtype = numpy.float64
    output_tdat_chunksize = 10000    # HDF5 chunksize for transition data (~500 KiB)
    tdat_buffersize = 1000000         # Internal buffer length (~50 MiB)
    
    def __init__(self, n_bins, output_group):
        self.n_bins = n_bins
        self.iibins = numpy.arange(n_bins)
        self.iibdisc = numpy.empty((n_bins,), numpy.bool_)
        
        self.bin_index_dtype = numpy.min_scalar_type(n_bins)
        
        self.tdat_dtype = numpy.dtype( [('timepoint',       self.index_dtype),
                                        ('initial_bin',     self.bin_index_dtype),
                                        ('final_bin',       self.bin_index_dtype),
                                        ('weight',          self.weight_dtype),
                                        ('initial_bin_pop', self.weight_dtype),
                                        ('duration',        self.index_dtype),
                                        ('fpt',             self.index_dtype),
                                        ])
        
        
        # HDF5 group in which to store results
        self.output_group = output_group
        self.tdat_buffer = numpy.empty((self.tdat_buffersize,), dtype=self.tdat_dtype)
        self.tdat_buffer_offset = 0
        self.output_tdat_offset = 0
        self.output_tdat_ds = None
        self.record_all_crossings = False
        self.record_self_transitions = False
                        
        # Accumulators/counters
        self.n_trans           = None # shape (n_bins,n_bins)

        # Time points and per-timepoint data
        self.last_exit         = None # (n_bins,)
        self.last_entry        = None # (n_bins,)
        self.last_completion   = None # (n_bins,n_bins) 
        self.bin_pops_last_exit = None # (n_bins,)       
        
        # Analysis continuation information
        self.timepoint          = None # current time index for separate calls on same trajectory
        self.last_bin           = None # last region occupied, for separate calls on same trajectory
        self.last_bin_pop       = None # total weight in self.last_region at end of last processing step
                
        self.clear()
        
    def clear(self):
        self.clear_state()
        self.n_trans         = numpy.zeros((self.n_bins,self.n_bins), self.count_dtype)
        self.tdat_buffer = numpy.empty((self.tdat_buffersize,), dtype=self.tdat_dtype)
        self.tdat_buffer_offset = 0
        self.output_tdat_offset = 0
        self.output_tdat_ds = None
        
    def clear_state(self):
        self.last_exit          = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_entry         = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_completion    = numpy.zeros((self.n_bins,self.n_bins), self.index_dtype)
        self.bin_pops_last_exit = numpy.zeros((self.n_bins,), self.weight_dtype)        
        self.timepoint          = 0
        self.last_bin           = None
        self.last_bin_pop       = None        

    def get_state(self):
        return {'last_entry':           self.last_entry.copy(),
                'last_exit':            self.last_exit.copy(),
                'last_completion':      self.last_completion.copy(),
                'bin_pops_last_exit':   self.bin_pops_last_exit.copy(),
                'timepoint':            self.timepoint,
                'last_bin':             self.last_bin,
                'last_bin_pop':         self.last_bin_pop,
                }
        
    def set_state(self, state_dict):
        self.last_entry = state_dict['last_entry']
        self.last_exit = state_dict['last_exit']
        self.last_completion = state_dict['last_completion']
        self.bin_pops_last_exit = state_dict['bin_pops_last_exit']
        self.timepoint = state_dict['timepoint']
        self.last_bin = state_dict['last_bin']
        self.last_bin_pop = state_dict['last_bin_pop']
        
    def record_transition_data(self, tdat):
        """Update running statistics and write transition data to HDF5 (with buffering)"""        
        
        # Write out accumulated transition data
        if self.output_tdat_ds is None:        
            # Create dataset
            try:
                del self.output_group['transitions']
            except KeyError:
                pass
            
            self.output_tdat_ds = self.output_group.create_dataset('transitions', shape=(1,), 
                                                                   dtype=self.tdat_dtype, maxshape=(None,), 
                                                                   chunks=(self.output_tdat_chunksize,),)
            
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
        
    def start_accumulation(self, assignments, weights, bin_pops):
        self.clear_state()
        timepoints = numpy.arange(len(assignments))
        self._accumulate_transitions(timepoints, assignments, weights, bin_pops)
    
    def continue_accumulation(self, assignments, weights, bin_pops):
        aug_assign = numpy.empty((len(assignments)+1,), assignments.dtype)
        aug_assign[0] = self.last_bin
        aug_assign[1:] = assignments
        
        aug_weights = numpy.empty((len(weights)+1,), self.weight_dtype)
        aug_weights[0] = 0
        aug_weights[1:] = weights
        
        aug_pops = numpy.empty((len(bin_pops)+1, len(bin_pops[0])), self.weight_dtype)
        aug_pops[0,:] = 0
        aug_pops[0, self.last_bin] = self.last_bin_pop
        
        timepoints = numpy.arange(self.timepoint, self.timepoint+len(aug_assign))
        
        self._accumulate_transitions(timepoints, aug_assign, aug_weights, aug_pops)
        
    
    def _accumulate_transitions(self, timepoints, assignments, weights, bin_pops):
        tdat = []
        
        trans_occur = assignments[1:] != assignments[:-1]
        trans_ibin  = assignments[:-1][trans_occur]
        trans_fbin  = assignments[1:][trans_occur]
        trans_timepoints = timepoints[1:][trans_occur]
        trans_weights   = weights[1:][trans_occur] # arrival weights
        trans_ibinpops = bin_pops[:-1][trans_occur]
        
        last_exit = self.last_exit
        last_entry = self.last_entry
        last_completion = self.last_completion
        bin_pops_last_exit = self.bin_pops_last_exit
        n_trans = self.n_trans
        iibdisc = self.iibdisc
        iibins = self.iibins
        tdat_maxlen = self.tdat_buffersize / 10
        record_all_crossings = self.record_all_crossings
        record_self_transitions = self.record_self_transitions
        for (trans_ti, weight, ibin, fbin, ibinpops) in izip(trans_timepoints, trans_weights, trans_ibin, trans_fbin, trans_ibinpops):
            # Record this crossing event's data
            
            if record_all_crossings:
                tdat.append((trans_ti, ibin, fbin, weight, ibinpops[ibin], 0, 0))
            
            iibdisc[:] = last_exit > 0
            iibdisc[ibin] = False
            iibdisc &= last_entry > last_completion[:,fbin]
            
            for iibin in iibins[iibdisc]:                
                duration = trans_ti - last_exit[iibin]
                if last_completion[iibin,fbin] > 0:
                    fpt      = trans_ti - last_completion[iibin,fbin]
                else:
                    fpt = 0

                if iibin != fbin or record_self_transitions:
                    tdat.append((trans_ti, iibin, fbin, weight, bin_pops_last_exit[iibin], duration, fpt))
                last_completion[iibin,fbin] = trans_ti
                n_trans[iibin,fbin] += 1

            last_exit[ibin] = trans_ti
            last_entry[fbin] = trans_ti
            last_completion[ibin,fbin] = trans_ti
            n_trans[ibin,fbin] += 1
            bin_pops_last_exit[ibin] = ibinpops[ibin]
            
            if len(tdat) > tdat_maxlen:
                self.record_transition_data(tdat)
        
        self.record_transition_data(tdat)
        self.timepoint = timepoints[-1]
        self.last_bin = assignments[-1]
        self.last_bin_pop = bin_pops[-1,assignments[-1]]
        