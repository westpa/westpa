# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.


import logging

log = logging.getLogger(__name__)

import numpy


import westpa
from oldtools.aframe import AnalysisMixin
from oldtools.aframe.trajwalker import TrajWalker

class TransitionEventAccumulator:
    index_dtype  = numpy.uintp
    count_dtype  = numpy.uint64
    weight_dtype = numpy.float64
    output_tdat_chunksize = 4096      # HDF5 chunksize for transition data (~300 KiB)
    tdat_buffersize = 524288          # Internal buffer length (~38 MiB)
    max_acc = 32768
    
    def __init__(self, n_bins, output_group, calc_fpts = True):
        self.calc_fpts = calc_fpts
        self.n_bins = n_bins
        self.iibins = numpy.arange(n_bins)
        self.iibdisc = numpy.empty((n_bins,), numpy.bool_)
        
        self.bin_index_dtype = numpy.min_scalar_type(n_bins)
        
        self.tdat_dtype = numpy.dtype( [('traj',             self.index_dtype),
                                        ('n_iter',           self.index_dtype),
                                        ('timepoint',        self.index_dtype),
                                        ('initial_bin',      self.bin_index_dtype),
                                        ('final_bin',        self.bin_index_dtype),
                                        ('initial_weight',   self.weight_dtype),
                                        ('final_weight',     self.weight_dtype),
                                        ('initial_bin_pop',  self.weight_dtype),
                                        ('duration',         self.index_dtype),
                                        ('fpt',              self.index_dtype),
                                        ])
        
        
        # HDF5 group in which to store results
        self.output_group = output_group
        self.tdat_buffer = numpy.empty((self.tdat_buffersize,), dtype=self.tdat_dtype)
        self.tdat_buffer_offset = 0
        self.output_tdat_offset = 0
        self.output_tdat_ds = None
                        
        # Accumulators/counters
        self.n_trans           = None # shape (n_bins,n_bins)

        # Time points and per-timepoint data
        self.last_exit          = None # (n_bins,)
        self.last_entry         = None # (n_bins,)
        self.last_completion    = None # (n_bins,n_bins)
        self.weight_last_exit   = None # (n_bins) 
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
        self.weight_last_exit   = numpy.zeros((self.n_bins,), self.weight_dtype)
        self.bin_pops_last_exit = numpy.zeros((self.n_bins,), self.weight_dtype)        
        self.timepoint          = 0
        self.last_bin           = None
        self.last_bin_pop       = None        

    def get_state(self):
        return {'last_entry':           self.last_entry.copy(),
                'last_exit':            self.last_exit.copy(),
                'last_completion':      self.last_completion.copy(),
                'weight_last_exit':     self.weight_last_exit.copy(),
                'bin_pops_last_exit':   self.bin_pops_last_exit.copy(),
                'timepoint':            self.timepoint,
                'last_bin':             self.last_bin,
                'last_bin_pop':         self.last_bin_pop,
                }
        
    def set_state(self, state_dict):
        self.last_entry = state_dict['last_entry']
        self.last_exit = state_dict['last_exit']
        self.last_completion = state_dict['last_completion']
        self.weight_last_exit = state_dict['weight_last_exit']
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
                                                                   chunks=(self.output_tdat_chunksize,),
                                                                   compression='gzip')
            
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
        
    def start_accumulation(self, assignments, weights, bin_pops, traj=0, n_iter=0):
        self.clear_state()
        timepoints = numpy.arange(len(assignments))
        self._accumulate_transitions(timepoints, assignments, weights, bin_pops, traj, n_iter)
    
    def continue_accumulation(self, assignments, weights, bin_pops, traj=0, n_iter=0):
        aug_assign = numpy.empty((len(assignments)+1,), assignments.dtype)
        aug_assign[0] = self.last_bin
        aug_assign[1:] = assignments
        
        aug_weights = numpy.empty((len(weights)+1,), self.weight_dtype)
        aug_weights[0] = 0
        aug_weights[1:] = weights
        
        aug_pops = numpy.empty((len(bin_pops)+1, len(bin_pops[0])), self.weight_dtype)
        aug_pops[0,:] = 0
        aug_pops[0, self.last_bin] = self.last_bin_pop
        aug_pops[1:] = bin_pops
        
        timepoints = numpy.arange(self.timepoint, self.timepoint+len(aug_assign))        
        self._accumulate_transitions(timepoints, aug_assign, aug_weights, aug_pops, traj, n_iter)
    
    def _accumulate_transitions(self, timepoints, assignments, weights, bin_pops, traj, n_iter):        
        tdat = []
        
        assignments_from_1 = assignments[1:]
        assignments_to_1   = assignments[:-1]
        
        calc_fpts = self.calc_fpts
        
        trans_occur = assignments_from_1 != assignments_to_1
        trans_ibin  = assignments_to_1[trans_occur]
        trans_fbin  = assignments_from_1[trans_occur]
        trans_timepoints = timepoints[1:][trans_occur]
        trans_weights   = weights[1:][trans_occur] # arrival weights
        trans_ibinpops = bin_pops[:-1][trans_occur]
        
        last_exit = self.last_exit
        last_entry = self.last_entry
        last_completion = self.last_completion
        bin_pops_last_exit = self.bin_pops_last_exit
        weight_last_exit = self.weight_last_exit
        n_trans = self.n_trans
        iibdisc = self.iibdisc
        iibins = self.iibins
        tdat_maxlen = self.max_acc
        for (trans_ti, weight, ibin, fbin, ibinpops) in zip(trans_timepoints, trans_weights, 
                                                             trans_ibin, trans_fbin, trans_ibinpops):
            # Record this crossing event's data
            bin_pops_last_exit[ibin] = ibinpops[ibin]            
            last_exit[ibin] = trans_ti
            last_entry[fbin] = trans_ti
            weight_last_exit[ibin] = weight
            
            # See what other transitions this crossing event completes
            iibdisc[:] = last_exit > 0
            iibdisc &= last_entry > last_completion[:,fbin]
            
            # Calculate event durations, etc for each transition generated by this crossing event
            durations = -last_exit + trans_ti + 1 # = time now - time of exit from initial bin
            if calc_fpts:
                fpts      = -last_completion[fbin,:] + trans_ti # = time now - time of last final->initial transition
                fpts[last_completion[:,fbin]==0] = 0
            
                for iibin in iibins[iibdisc]:
                    tdat.append((traj, n_iter,trans_ti,iibin,fbin,weight_last_exit[iibin], 
                                 weight, bin_pops_last_exit[iibin],durations[iibin], fpts[iibin]))
            else:
                for iibin in iibins[iibdisc]:
                    tdat.append((traj, n_iter,trans_ti,iibin,fbin,weight_last_exit[iibin], 
                                 weight, bin_pops_last_exit[iibin],durations[iibin], 0))    
            
            # Update tracking and statistics                    
            last_completion[iibdisc,fbin] = trans_ti
            n_trans[iibdisc,fbin] += 1
            
            if len(tdat) > tdat_maxlen:
                self.record_transition_data(tdat)
        
        self.record_transition_data(tdat)
        self.timepoint = timepoints[-1]
        self.last_bin = assignments[-1]
        self.last_bin_pop = bin_pops[-1,assignments[-1]]

class TransitionAnalysisMixin(AnalysisMixin):
    def __init__(self):
        super(TransitionAnalysisMixin,self).__init__()
        self.discard_transition_data = False
        self.calc_fpts = False
        self.trans_h5gname = 'transitions'
        self.trans_h5group = None
        self.__transitions_ds = None
        self.n_trajs = 0
        
    def require_transitions_group(self):
        if self.trans_h5group is None:
            self.trans_h5group = self.anal_h5file.require_group(self.trans_h5gname)
        return self.trans_h5group

    def delete_transitions_group(self):
        self.trans_h5group = None
        del self.anal_h5file[self.trans_h5gname]
                
    def get_transitions_ds(self):
        if self.__transitions_ds is not None:
            return self.__transitions_ds
        else:
            self.__transitions_ds = self.trans_h5group['transitions']
            return self.__transitions_ds
    
    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc = super(TransitionAnalysisMixin,self).add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
        
        group = parser.add_argument_group('transition analysis options')
        group.add_argument('--discard-transition-data', dest='discard_transition_data', action='store_true',
                           help='''Discard any existing transition data stored in the analysis HDF5 file.''')        
    
    def process_args(self, args, upcall = True):                
        self.discard_transition_data = args.discard_transition_data
        
        
        if upcall:
            try:
                upfunc = super(TransitionAnalysisMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)
                
    def require_transitions(self):
        self.require_bin_assignments()
        self.require_transitions_group()
        do_trans = False
        if self.discard_transition_data:
            westpa.rc.pstatus('Discarding existing transition data.')
            do_trans = True
        elif not self.check_data_binhash(self.trans_h5group):
            westpa.rc.pstatus('Bin definitions have changed; deleting existing transition data.')
            do_trans = True
        elif 'transitions' in self.trans_h5group and not self.check_data_iter_range_least(self.trans_h5group): 
            westpa.rc.pstatus('Existing transition data is for different first/last iterations; deleting.')
            do_trans = True
                    
        if do_trans:
            self.delete_transitions_group()
            self.find_transitions()

    def find_transitions(self):
        westpa.rc.pstatus('Finding transitions...')
        output_group = self.require_transitions_group()
            
        self.n_segs_visited = 0
        self.n_total_segs = self.total_segs_in_range(self.first_iter,self.last_iter)
        self.accumulator = TransitionEventAccumulator(self.n_bins, output_group, calc_fpts = self.calc_fpts)        
        self.bin_assignments = self.get_bin_assignments(self.first_iter,self.last_iter)
        self.bin_populations = self.get_bin_populations(self.first_iter,self.last_iter)
        
        walker = TrajWalker(data_reader = self)
        
        self.__pcoord_len = self.get_pcoord_len(self.first_iter)
        self.__quiet_mode = westpa.rc.quiet_mode
        
        walker.trace_trajectories(self.first_iter, self.last_iter, callable=self._segment_callback, include_pcoords=False,
                                  get_state = self.accumulator.get_state, set_state = self.accumulator.set_state)
        self.accumulator.flush_transition_data()
        try:
            del output_group['n_trans']
        except KeyError:
            pass
        output_group['n_trans'] = self.accumulator.n_trans
        
        for h5object in (output_group, output_group['n_trans'], output_group['transitions']):
            self.record_data_iter_range(h5object)
            self.record_data_binhash(h5object)
            h5object.attrs['n_trajs'] = self.n_trajs
            h5object.attrs['n_segs'] = self.n_segs_visited
        
        self.accumulator.clear()
        westpa.rc.pstatus()
        
    def _segment_callback(self, segment, children, history):
        iiter = segment.n_iter - self.first_iter
        seg_id = segment.seg_id
        weights = numpy.empty((self.__pcoord_len,), numpy.float64)
        weights[:] = segment.weight
        bin_pops = self.bin_populations[iiter, :, :]
        
        if len(history) == 0:
            # New trajectory
            self.n_trajs += 1
            self.accumulator.start_accumulation(self.bin_assignments[iiter, seg_id, :], weights, bin_pops, 
                                                traj=self.n_trajs, n_iter=segment.n_iter)
        else:
            # Continuing trajectory
            self.accumulator.continue_accumulation(self.bin_assignments[iiter, seg_id, :], weights, bin_pops, 
                                                   traj=self.n_trajs, n_iter=segment.n_iter)
            
        self.n_segs_visited += 1
        
        if not self.__quiet_mode and (self.n_segs_visited % 1000 == 0 or self.n_segs_visited == self.n_total_segs):
            pct_visited = self.n_segs_visited / self.n_total_segs * 100
            westpa.rc.pstatus('\r  {:d} of {:d} segments ({:.1f}%) analyzed ({:d} independent trajectories)'
                            .format(int(self.n_segs_visited), int(self.n_total_segs), float(pct_visited), self.n_trajs), 
                            end='')
            westpa.rc.pflush()

class BFTransitionAnalysisMixin(TransitionAnalysisMixin):
    
    def require_transitions(self):
        self.require_bin_assignments()
        self.require_transitions_group()
        do_trans = False
        if self.discard_transition_data:
            westpa.rc.pstatus('Discarding existing transition data.')
            do_trans = True
        elif not self.check_data_binhash(self.trans_h5group):
            westpa.rc.pstatus('Bin definitions have changed; deleting existing transition data.')
            do_trans = True
        
        if do_trans:
            self.delete_transitions_group()
            self.find_transitions()

    def find_transitions(self,chunksize=65536):
        self.require_bf_h5file()
        self.require_binning_group()
        westpa.rc.pstatus('Finding transitions...')
        output_group = self.require_analysis_group('transitions')
            
        self.accumulator = TransitionEventAccumulator(self.n_bins, output_group, calc_fpts = True)
        assignments_ds = self.binning_h5group['bin_assignments']
        
        max_nrows = assignments_ds.len()
        maxwidth_nrows = len(str(max_nrows))
        
        for traj_id in range(self.get_n_trajs()):
            nrows = self.get_traj_len(traj_id)
            for istart in range(0, nrows, chunksize):
                iend = min(istart+chunksize,nrows)
                assignments = assignments_ds[traj_id,istart:iend]
                weights = numpy.ones((len(assignments),))
                binpops = numpy.ones((len(assignments),self.n_bins))
                
                if istart == 0:
                    self.accumulator.start_accumulation(assignments, weights, binpops, traj=traj_id)
                else:
                    self.accumulator.continue_accumulation(assignments, weights, binpops, traj=traj_id)
                    
                westpa.rc.pstatus('\r  Trajectory {:d}: {:{mwnr}d}/{:<{mwnr}d} ({:.2f}%)'
                                .format(int(traj_id), int(iend), int(nrows), iend/nrows*100,mwnr=maxwidth_nrows), 
                                end='')
                westpa.rc.pflush()
                self.accumulator.flush_transition_data()
                del assignments, weights, binpops
            westpa.rc.pstatus()                 
                
        try:
            del output_group['n_trans']
        except KeyError:
            pass
        output_group['n_trans'] = self.accumulator.n_trans
        
        for h5object in (output_group, output_group['n_trans'], output_group['transitions']):
            self.record_data_binhash(h5object)
        
        self.accumulator.clear()
