from __future__ import division; __metaclass__ = type

import sys, time, operator, math, numpy, re, random
from itertools import izip
from datetime import timedelta
from collections import namedtuple
import logging
log = logging.getLogger(__name__)

import wemd
from wemd.states import BasisState, InitialState, TargetState
from wemd.util import extloader
from wemd import Segment
from wemd.util.miscfn import vgetattr

from wemd.work_managers import ops as wm_ops
from wemd.data_manager import weight_dtype

EPS = numpy.finfo(numpy.float64).eps

RecyclingInfo = namedtuple('RecyclingInfo', ['count', 'weight'])

class PropagationError(RuntimeError):
    pass 

class WESimManager:
    def __init__(self):        
        self.data_manager = wemd.rc.data_manager
        self.we_driver = wemd.rc.we_driver
        self.work_manager = wemd.rc.work_manager
        self.propagator = wemd.rc.propagator
        self.system = wemd.rc.system
                                
        # A table of function -> list of (priority, name, callback) tuples
        self._callback_table = {}
        self._valid_callbacks = set((self.prepare_run, self.finalize_run,
                                     self.prepare_iteration, self.finalize_iteration,
                                     self.pre_propagation, self.post_propagation,
                                     self.pre_we, self.post_we))
        self._callbacks_by_name = {fn.__name__: fn for fn in self._valid_callbacks}
        self._n_propagated = 0
        
        self.do_gen_istates = wemd.rc.config.get_bool('system.gen_istates', False) 
        
        # Per-iteration variables
        self.n_iter = None                  # current iteration
        self.target_states = None           # TargetStates valid at this iteration
        self.target_state_bins = None       # Bins (in self.final_binning) associated with each TargetState
        
        self.current_iter_bstates = None       # BasisStates valid at this iteration
        self.current_iter_istates = None       # InitialStates used in this iteration        
        self.next_iter_bstates = None          # BasisStates valid for the next iteration
        self.next_iter_bstate_cprobs = None    # Cumulative probabilities for basis states, used for selection
        self.next_iter_istates = None
        self.next_iter_spare_istates = None    # InitialStates available for use next iteration
        self.next_iter_assigned_istates = None # InitialStates that were available or generated in this iteration but then used
        
        self.initial_binning = None         # Binning of segments at beginning of this iteration
        self.final_binning = None           # Binning of segments at the end of this iteration
        self.segments = None                # Mapping of seg_id to segment for all segments in this iteration
        self.completed_segments = None      # Mapping of seg_id to segment for all completed segments in this iteration
        self.incomplete_segments = None     # Mapping of seg_id to segment for all incomplete segments in this iteration
        self.to_recycle = None              # Mapping of seg_id to segment for all completed segments to be recycled
        self.n_bins = None
        
        self.next_iter_segments = None      # List of segments to be created for the next iteration
        self.next_iter_binning = None       # Binning of next iteration's segments
        
    def register_callback(self, hook, function, priority=0):
        '''Registers a callback to execute during the given ``hook`` into the simulation loop. The optional
        priority is used to order when the function is called relative to other registered callbacks.'''
        
        if hook not in self._valid_callbacks:
            try:
                hook = self._callbacks_by_name[hook]
            except KeyError:
                raise KeyError('invalid hook {!r}'.format(hook))
            
        try:
            self._callback_table[hook].add((priority,function.__name__,function))
        except KeyError:
            self._callback_table[hook] = set([(priority,function.__name__,function)])
        
        log.debug('registered callback {!r} for hook {!r}'.format(function, hook))
                
    def invoke_callbacks(self, hook, *args, **kwargs):
        callbacks = self._callback_table.get(hook, [])
        sorted_callbacks = sorted(callbacks)
        for (priority, name, fn) in sorted_callbacks:
            fn(*args, **kwargs)
    
    def load_plugins(self):
        plugin_text = wemd.rc.config.get('plugins.enable','')
        plugin_names = re.split(r'\s*,\s*', plugin_text.strip())
        for plugin_name in plugin_names:
            if not plugin_name: continue
            
            log.info('loading plugin {!r}'.format(plugin_name))
            plugin = extloader.get_object(plugin_name)(self)
            log.debug('loaded plugin {!r}'.format(plugin))

    def report_bin_statistics(self, region_set, save_summary=False):
        segments = self.segments.values()
        bins = region_set.get_all_bins()
        bin_counts = vgetattr('count', bins, numpy.uint)
        target_counts = vgetattr('target_count', bins, numpy.uint)

        # Do not include bins with target count zero (e.g. sinks, never-filled bins) in the (non)empty bins statistics
        n_active_bins = len(target_counts[target_counts!=0])
        seg_probs  = vgetattr('weight', segments, numpy.float64)
        bin_probs  = vgetattr('weight', bins, numpy.float64)
        norm = seg_probs.sum()
        
        assert abs(1 - norm) < EPS*(len(segments)+n_active_bins)
        
        min_seg_prob = seg_probs[seg_probs!=0].min()
        max_seg_prob = seg_probs.max()
        seg_drange   = math.log(max_seg_prob/min_seg_prob)
        min_bin_prob = bin_probs[bin_probs!=0].min()
        max_bin_prob = bin_probs.max()
        bin_drange = math.log(max_bin_prob/min_bin_prob)
        n_pop = len(bin_counts[bin_counts!=0])
        
        wemd.rc.pstatus('{:d} of {:d} ({:%}) active bins are populated'.format(n_pop, n_active_bins,n_pop/n_active_bins))
        wemd.rc.pstatus('per-bin minimum non-zero probability:       {:g}'.format(min_bin_prob))
        wemd.rc.pstatus('per-bin maximum probability:                {:g}'.format(max_bin_prob))
        wemd.rc.pstatus('per-bin probability dynamic range (kT):     {:g}'.format(bin_drange))
        wemd.rc.pstatus('per-segment minimum non-zero probability:   {:g}'.format(min_seg_prob))
        wemd.rc.pstatus('per-segment maximum non-zero probability:   {:g}'.format(max_seg_prob))
        wemd.rc.pstatus('per-segment probability dynamic range (kT): {:g}'.format(seg_drange))
        wemd.rc.pstatus('norm = {:g}, error in norm = {:g} ({:.2g}*epsilon)'.format(norm,(norm-1),(norm-1)/EPS))
        wemd.rc.pflush()
        
        if save_summary:
            iter_summary = self.data_manager.get_iter_summary()
            iter_summary['n_particles'] = len(segments)
            iter_summary['norm'] = norm
            iter_summary['min_bin_prob'] = min_bin_prob
            iter_summary['max_bin_prob'] = max_bin_prob
            iter_summary['min_seg_prob'] = min_seg_prob
            iter_summary['max_seg_prob'] = max_seg_prob
            if numpy.isnan(iter_summary['cputime']): iter_summary['cputime'] = 0.0
            if numpy.isnan(iter_summary['walltime']): iter_summary['walltime'] = 0.0
            self.data_manager.update_iter_summary(iter_summary)

    def get_bstate_pcoords(self, basis_states):
        '''For each of the given ``basis_states``, calculate progress coordinate values
        as necessary.  The HDF5 file is not updated.'''
        
        wemd.rc.pstatus('Calculating progress coordinate values for basis states.')
        futures = [self.work_manager.submit(wm_ops.get_pcoord, self.propagator, basis_state)
                   for basis_state in basis_states]
        fmap = {future: i for (i, future) in enumerate(futures)}
        for future in self.work_manager.as_completed(futures): 
            basis_states[fmap[future]].pcoord = future.get_result().pcoord
        
    def report_basis_states(self, basis_states):
        pstatus = wemd.rc.pstatus
        pstatus('{:d} basis state(s) present'.format(len(basis_states)), end='')
        if wemd.rc.verbose_mode:
            pstatus(':')
            pstatus('{:6s}    {:12s}    {:20s}    {:20s}    {}'
                             .format('ID', 'Label', 'Probability', 'Aux Reference', 'Progress Coordinate'))
            for basis_state in basis_states:
                pstatus('{:<6d}    {:12s}    {:<20.14g}    {:20s}    {}'.
                                 format(basis_state.state_id, basis_state.label, basis_state.probability, basis_state.auxref or '',
                                        ', '.join(map(str,basis_state.pcoord))))
        pstatus()
        wemd.rc.pflush()
        
    def report_target_states(self, target_states):
        pstatus = wemd.rc.pstatus
        pstatus('{:d} target state(s) present'.format(len(target_states)), end='')    
        if wemd.rc.verbose_mode and target_states:
            pstatus(':')
            pstatus('{:6s}    {:12s}    {}'.format('ID', 'Label', 'Progress Coordinate'))
            for target_state in target_states:
                pstatus('{:<6d}    {:12s}    {}'
                                 .format(target_state.state_id, target_state.label, ','.join(map(str,target_state.pcoord))))        
        pstatus()
        wemd.rc.pflush()
        
    def initialize_simulation(self, basis_states, target_states, segs_per_state=1, suppress_we=False):
        '''Initialize a new weighted ensemble simulation, taking ``segs_per_state`` initial
        states from each of the given ``basis_states``.
        
        ``w_init`` is the forward-facing version of this function'''
        
        data_manager = self.data_manager
        work_manager = self.work_manager
        propagator = self.propagator
        pstatus = wemd.rc.pstatus
        system = self.system

        pstatus('Creating HDF5 file {!r}'.format(self.data_manager.we_h5filename))
        data_manager.prepare_backing()

        # Process target states
        data_manager.save_target_states(target_states)
        self.report_target_states(target_states)

        
        # Process basis states
        self.get_bstate_pcoords(basis_states)        
        self.data_manager.create_ibstate_group(basis_states)
        self.report_basis_states(basis_states)
        
        pstatus('Preparing initial states')
        segments = []
        initial_states = []
        if self.do_gen_istates:
            istate_type = InitialState.ISTATE_TYPE_GENERATED
        else:
            istate_type = InitialState.ISTATE_TYPE_BASIS
            
        for basis_state in basis_states:
            for iseg in xrange(segs_per_state):
                initial_state = data_manager.create_initial_states(1,1)[0]
                initial_state.basis_state_id =  basis_state.state_id
                initial_state.basis_state = basis_state
                initial_state.istate_type = istate_type
                segment = Segment(weight=basis_state.probability/segs_per_state,pcoord=system.new_pcoord_array(),
                                  parent_id=-(initial_state.state_id+1), wtg_parent_ids=(-(initial_state.state_id+1),),
                                  )
                initial_states.append(initial_state)
                log.debug('initial state created: {!r}'.format(initial_state))
                segments.append(segment)
                
        if self.do_gen_istates:
            futures = [work_manager.submit(wm_ops.gen_istate, propagator, initial_state.basis_state, initial_state)
                       for initial_state in initial_states]
            for future in work_manager.as_completed(futures):
                rbstate, ristate = future.get_result()
                initial_states[ristate.state_id].pcoord = ristate.pcoord
                segments[ristate.state_id].pcoord[0] = ristate.pcoord
        else:
            for segment, initial_state in izip(segments, initial_states):
                initial_state.pcoord = basis_state.pcoord
                initial_state.istate_status = InitialState.ISTATE_STATUS_PREPARED
                segment.pcoord[0] = basis_state.pcoord

        tprob = sum(segment.weight for segment in segments)
        if abs(1.0 - tprob) > len(segments) * EPS:
            pscale = 1.0/tprob
            log.warning('Weights of initial segments do not sum to unity; scaling by {:g}'.format(pscale))
            for segment in segments:
                segment.weight *= pscale
                    
        data_manager.update_initial_states(initial_states, n_iter=1)
        
        if not suppress_we:
            # TODO: what to do if something winds up in a recycling region?
            # At the moment, fail with an exception
            region_set = system.new_region_set()
            new_region_set = self.we_driver.run_we(region_set, segments)
        else:
            new_region_set = system.new_region_set()
            new_region_set.assign_to_bins(segments, key=Segment.initial_pcoord)
        
        all_bins = new_region_set.get_all_bins()

        if target_states:
            for target_index in new_region_set.map_to_all_indices([target.pcoord for target in target_states]):
                all_bins[target_index].target_count = 0

        bin_occupancies = numpy.array(map(operator.attrgetter('count'), all_bins))
        target_occupancies = numpy.array(map(operator.attrgetter('target_count'), all_bins))
        segments = list(new_region_set.particles)

        # Make sure we have a norm of 1
        for segment in segments:
            segment.n_iter = 1
            segment.status = Segment.SEG_STATUS_PREPARED
            assert segment.parent_id < 0
            initial_states[segment.initial_state_id].iter_used = 1        
                    
                    
        data_manager.prepare_iteration(1, segments)
        data_manager.update_initial_states(initial_states, n_iter=1)
                    
        if wemd.rc.verbose_mode:
            pstatus('\nSegments generated:')
            for segment in segments:
                pstatus('{!r}'.format(segment))
        
        
        pstatus('''
        Total bins:            {total_bins:d}
        Initial replicas:      {init_replicas:d} in {occ_bins:d} bins, total weight = {weight:g}
        Total target replicas: {total_replicas:d}
        '''.format(total_bins=len(all_bins), init_replicas=sum(bin_occupancies), occ_bins=len(bin_occupancies[bin_occupancies > 0]),
                   weight = sum(segment.weight for segment in segments), total_replicas = sum(target_occupancies)))
        
        # Send the segments over to the data manager to commit to disk            
        data_manager.current_iteration = 1
        
        # Report statistics
        pstatus('Simulation prepared.')
        self.segments = {segment.seg_id: segment for segment in segments}
        self.report_bin_statistics(new_region_set,save_summary=True)
        data_manager.flush_backing()

    def prepare_iteration(self):
        log.debug('beginning iteration {:d}'.format(self.n_iter))
        
        # Clean up from last iteration
        # Explicit deletes are used to make sure that two iterations' worth of segment data don't wind up in RAM at once
        del self.segments, self.next_iter_segments, self.initial_binning, self.final_binning, self.next_iter_binning
        self.next_iter_segments = None
        self.next_iter_binning = None
        
        # Prepare region sets (bins) for this iteration
        # Since we track transitions, it's easiest to have separate region sets for initial and final states of each segment
        self.initial_binning = self.system.new_region_set()
        self.final_binning = self.system.new_region_set()
        self.n_bins = len(self.initial_binning.get_all_bins())
                
        # Store bin identity hash in HDF5 to detect when bins have changed
        # We directly modify an HDF5 object (the 'binhash' attribute on the iteration group), so obtain the lock first
        binhash = self.initial_binning.identity_hash().hexdigest()
        with self.data_manager.lock:
            iter_summary = self.data_manager.get_iter_summary(self.n_iter)
            iter_summary['binhash'] = binhash
            iter_group = self.data_manager.get_iter_group(self.n_iter)
            iter_group.attrs['binhash'] = binhash
            self.data_manager.update_iter_summary(iter_summary, self.n_iter)
        
        # Get target states and map them to bins
        self.target_states = self.data_manager.get_target_states(self.n_iter)

        if self.target_states:
            self.target_state_bins = list(self.final_binning.map_to_bins([target_state.pcoord for target_state in self.target_states]))
            for bin in self.target_state_bins:
                bin.target_count = 0
            log.debug('target_state_bins={!r}'.format(self.target_state_bins))
        
        # Get basis states used in this iteration
        self.current_iter_bstates = self.data_manager.get_basis_states(self.n_iter)
        
        # Get the segments for this iteration and separate into complete and incomplete
        segments = self.segments = {segment.seg_id: segment for segment in self.data_manager.get_segments()}
        log.debug('loaded {:d} segments'.format(len(segments)))
        completed_segments = self.completed_segments = {}
        incomplete_segments = self.incomplete_segments = {}
        for segment in segments.itervalues():
            if segment.status == Segment.SEG_STATUS_COMPLETE:
                completed_segments[segment.seg_id] = segment
            else:
                incomplete_segments[segment.seg_id] = segment
        log.debug('{:d} segments are complete; {:d} are incomplete'.format(len(completed_segments), len(incomplete_segments)))
        
        if len(incomplete_segments) == len(segments):
            # Starting a new iteration
            wemd.rc.pstatus('Beginning iteration {:d}'.format(self.n_iter))
        elif incomplete_segments:
            wemd.rc.pstatus('Continuing iteration {:d}'.format(self.n_iter))
        wemd.rc.pstatus('{:d} segments remain in iteration {:d} ({:d} total)'.format(len(incomplete_segments), self.n_iter,
                                                                                     len(segments)))
        
        # Get the initial states active for this iteration (so that the propagator has them if necessary)
        self.current_iter_istates = {state.state_id: state for state in 
                                     self.data_manager.get_segment_initial_states(segments.values())}
        log.debug('This iteration uses {:d} initial states'.format(len(self.current_iter_istates)))
        
        # Assign this iteration's segments' initial points to bins
        self.initial_binning.assign_to_bins(segments.itervalues(), key=Segment.initial_pcoord)
        self.report_bin_statistics(self.initial_binning, save_summary=True)
        
        # Do the same for final bins, while also noting segments that wind up in recycling regions 
        self.to_recycle = {}
        if completed_segments:
            self.final_binning.assign_to_bins(completed_segments.values(), key=Segment.final_pcoord)
            
            if self.target_states:
                for bin, target_state in izip(self.target_state_bins, self.target_states):
                    log.debug('bin={!r}, target_state={!r}'.format(bin,target_state))
                    if len(bin):
                        log.debug('{:d} replicas in target state {!r}'.format(len(bin), target_state))
                        self.to_recycle.update({segment.seg_id: segment for segment in bin})
        log.debug('{:d} replicas in target states'.format(len(self.to_recycle)))
        
        # Get the basis states and initial states for the next iteration, necessary for doing on-the-fly recycling 
        self.next_iter_bstates = self.data_manager.get_basis_states(self.n_iter+1)
        self.next_iter_bstate_cprobs = numpy.add.accumulate([bstate.probability for bstate in self.next_iter_bstates])
        self.next_iter_assigned_istates = set()
        self.next_iter_spare_istates = set(self.data_manager.get_unused_initial_states(n_iter=self.n_iter+1))
        # No segments can exist for the next iteration yet, so this suffices to catch all valid states for the next iteration
        self.next_iter_istates = {state.state_id: state for state in self.next_iter_spare_istates}
        log.debug('{:d} unused initial states found'.format(len(self.next_iter_spare_istates)))
        
        # Invoke callbacks
        self.invoke_callbacks(self.prepare_iteration)
        
        log.debug('dispatching propagator prep_iter to work manager')        
        self.work_manager.submit(wm_ops.prep_iter, self.propagator, self.n_iter, segments).get_result()
        
    def finalize_iteration(self):
        '''Perform customized processing/cleanup on just-completed segments at the end of an iteration'''
        log.debug('finalizing iteration {:d}'.format(self.n_iter))
        
        self.invoke_callbacks(self.finalize_iteration)
        
        log.debug('dispatching propagator post_iter to work manager')
        self.work_manager.submit(wm_ops.post_iter, self.propagator, self.n_iter, self.segments.values()).get_result()


    def get_istate_futures(self, n_states):
        '''Add ``n_states`` initial states to the internal list of initial states assigned to
        recycled particles.  Spare states are used if available, otherwise new states are created.
        If created new initial states requires generation, then a set of futures is returned
        representing work manager tasks corresponding to the necessary generation work.'''
        
        futures = set()
        for i in xrange(n_states):
            if self.next_iter_spare_istates:
                log.debug('assigning spare state')
                self.next_iter_assigned_istates.add(self.next_iter_spare_istates.pop())
            else:
                # Select a basis state according to its weight
                ibstate = numpy.digitize([random.random()], self.next_iter_bstate_cprobs)
                basis_state = self.next_iter_bstates[ibstate]
                initial_state = self.data_manager.create_initial_states(1, n_iter=self.n_iter+1)[0]
                initial_state.iter_created = self.n_iter #TODO: this doesn't seem to fit with the above n_iter+1; make conformant?
                initial_state.basis_state_id = basis_state.state_id
                initial_state.istate_status = InitialState.ISTATE_STATUS_PENDING
                self.next_iter_istates[initial_state.state_id] = initial_state
                
                if self.do_gen_istates:
                    log.debug('generating new initial state from basis state {!r}'.format(basis_state))
                    initial_state.istate_type = InitialState.ISTATE_TYPE_GENERATED
                    futures.add(self.work_manager.submit(wm_ops.gen_istate,self.propagator, basis_state, initial_state))
                else:
                    log.debug('using basis state {!r} directly'.format(basis_state))
                    initial_state.istate_type = InitialState.ISTATE_TYPE_BASIS
                    initial_state.pcoord = basis_state.pcoord.copy()
                    self.assigned_initial_states.add(initial_state)
                self.data_manager.update_initial_states([initial_state], n_iter=self.n_iter+1)
        return futures
                                    
    def propagate(self):
        segments = self.incomplete_segments.values()
        log.debug('iteration {:d}: propagating {:d} segments'.format(self.n_iter, len(segments)))
        futures = set()        
        segment_futures = set()
        istate_gen_futures = self.get_istate_futures(len(self.to_recycle))
        futures.update(istate_gen_futures)
        
        log.debug('there are {:d} segments in target regions, which require generation of {:d} initial states'
                  .format(len(self.to_recycle),len(istate_gen_futures)))
                
        for segment in segments:
            pbstates, pistates = wemd.states.pare_basis_initial_states(self.current_iter_bstates, 
                                                                       self.current_iter_istates.values(), [segment])
            future = self.work_manager.submit(wm_ops.propagate, self.propagator, pbstates, pistates, [segment])
            futures.add(future)
            segment_futures.add(future)
        
        while futures:
            future = self.work_manager.wait_any(futures)
            futures.remove(future)
            
            if future in segment_futures:
                segment_futures.remove(future)
                incoming = future.get_result()
                self._n_propagated += 1
                log.debug('recording results for {!r}, {:d} for this run'.format(incoming, self._n_propagated))
                
                self.segments.update({segment.seg_id: segment for segment in incoming})
                self.completed_segments.update({segment.seg_id: segment for segment in incoming})
                
                for segment in incoming:
                    final_bin = self.final_binning.map_to_bins([segment.pcoord[-1]])[0]
                    final_bin.add(segment)
                    log.debug('incoming segment {!r} mapped to bin {!r}'.format(segment, final_bin))

                    for target_bin in self.target_state_bins or []:
                        if final_bin is target_bin:
                            # This particle must be recycled
                            log.debug('segment {!r} will be recycled (assigned to {!r})'.format(segment, final_bin))
                            self.to_recycle[segment.seg_id] = segment
                            new_futures = self.get_istate_futures(1)
                            if new_futures:
                                log.debug('new futures created for initial state generation: {!r}'.format(new_futures))
                            istate_gen_futures.update(new_futures)
                            futures.update(new_futures)
                            break

                with self.data_manager.flushing_lock():                        
                    self.data_manager.update_segments(self.n_iter, incoming)

            elif future in istate_gen_futures:
                istate_gen_futures.remove(future)
                _basis_state, initial_state = future.get_result()
                log.debug('received newly-prepared initial state {!r}'.format(initial_state))
                initial_state.istate_status = InitialState.ISTATE_STATUS_PREPARED
                self.data_manager.update_initial_states([initial_state], n_iter=self.n_iter+1)
                self.next_iter_assigned_istates.add(initial_state)
                self.next_iter_istates[initial_state.state_id] = initial_state
                self.data_manager.flush_backing()
            else:
                raise AssertionError('untracked future')                    
                    
        log.debug('done with propagation')
        self.save_bin_data()
        
    def save_bin_data(self):
        '''Calculate and write bin assignments, populations, transition counts, fluxes, etc to HDF5. The v0.7 code
        saved this on a per-timepoint basis, but this requires essentially two copies of the progress coordinate
        data for the entire set of segments to reside in RAM. This version saves only initial and final bin assignments,
        '''
        # save_bin_data(self, populations, n_trans, fluxes, rates, n_iter=None)
        
        n_segs = len(self.segments)
        n_bins = self.n_bins
        
        assignments = numpy.empty((n_segs, 2), numpy.min_scalar_type(n_bins))
        populations = numpy.zeros((n_bins,2), weight_dtype)
        n_trans = numpy.zeros((n_bins,n_bins), numpy.uint)
        fluxes = numpy.zeros((n_bins,n_bins), weight_dtype)
        rates = numpy.zeros((n_bins,n_bins), weight_dtype)
        
        # Though already assigned to bins in self.initial_binning and self.final_binning, it's
        # very likely faster to re-do the assignments than to search bins for segments
        seg_ids = sorted(self.segments.iterkeys())
        
        assignments[:,0] = self.initial_binning.map_to_all_indices([self.segments[seg_id].pcoord[0] for seg_id in seg_ids])
        assignments[:,-1] = self.final_binning.map_to_all_indices([self.segments[seg_id].pcoord[-1] for seg_id in seg_ids])
        
        segments = self.segments
        for seg_id, init_assignment, final_assignment in izip(seg_ids, assignments[:,0], assignments[:,-1]):
            segment = segments[seg_id]
            weight = segment.weight
            
            populations[init_assignment,0] += weight
            populations[final_assignment,-1] += weight
            n_trans[init_assignment,final_assignment] += 1
            fluxes[init_assignment,final_assignment] += weight
        
        for i in xrange(0,n_bins):
            if populations[i,0] > 0:
                rates[i,:] = fluxes[i,:] / populations[i,0] 
        
        self.data_manager.save_bin_data(assignments, populations, n_trans, fluxes, rates)
                

    def check_propagation(self):
        failed_segments = [segment for segment in self.segments.itervalues() if segment.status != Segment.SEG_STATUS_COMPLETE]
        
        if failed_segments:
            failed_ids = '  \n'.join(str(segment.seg_id) for segment in failed_segments)
            log.error('propagation failed for {:d} segment(s):\n{}'.format(len(failed_segments), failed_ids))
            raise PropagationError('propagation failed for {:d} segments'.format(len(failed_segments)))
        else:
            log.debug('propagation complete for iteration {:d}'.format(self.n_iter))
            
        failed_istates = [istate for istate in self.next_iter_assigned_istates
                          if istate.istate_status != InitialState.ISTATE_STATUS_PREPARED]
        log.debug('{!r}'.format(failed_istates))
        if failed_istates:
            failed_ids = '  \n'.join(str(istate.state_id) for istate in failed_istates)
            log.error('initial state generation failed for {:d} states:\n{}'.format(len(failed_istates), failed_ids))
            raise PropagationError('initial state generation failed for {:d} states'.format(len(failed_istates)))
        else:
            log.debug('initial state generation complete for iteration {:d}'.format(self.n_iter))


    def run_we(self):
        '''Run the weighted ensemble algorithm based on the binning in self.final_bins and
        the recycled particles in self.to_recycle, creating and committing the next iteration's
        segments to storage as well.'''

        # Remove recycled particles from target bins
        recycled_segs = set(self.to_recycle.values())
        for recycled_seg in recycled_segs:
            recycled_seg.endpoint_type = Segment.SEG_ENDPOINT_RECYCLED
        
        recycling_info = []

        if self.target_state_bins:
            for target_bin in self.target_state_bins:
                recycling_info.append(RecyclingInfo(len(target_bin), target_bin.weight))
                target_bin.difference_update(recycled_segs)
                assert len(target_bin) == 0
            self.data_manager.save_recycling_data(recycling_info)

        # assign initial states to particles being recycled
        init_segs = []
        used_initial_states = set()
        for segment, initial_state in izip(self.to_recycle.itervalues(), self.next_iter_assigned_istates):
            new_segment = Segment(n_iter=self.n_iter+1, parent_id=-(initial_state.state_id+1),
                                  weight=segment.weight, pcoord=self.system.new_pcoord_array(),
                                  wtg_parent_ids=[-(initial_state.state_id+1)],)
            new_segment.pcoord[0] = initial_state.pcoord
            init_segs.append(new_segment)
        
        log.debug('creating {:d} new particles due to recycling: {!r}'.format(len(init_segs), init_segs))                
        new_region_set = self.we_driver.run_we(self.final_binning, init_segs)
        new_segments = list(new_region_set.particles)
        
        for segment in new_segments:
            segment.n_iter = self.n_iter+1
            segment.status = Segment.SEG_STATUS_PREPARED
            if segment.initpoint_type == Segment.SEG_INITPOINT_NEWTRAJ:
                initial_state = self.next_iter_istates[segment.initial_state_id]
                initial_state.iter_used = self.n_iter+1
                used_initial_states.add(initial_state)            
                        
        if used_initial_states:
            self.data_manager.update_initial_states(used_initial_states)
            
        self.next_iter_segments = new_segments
        self.next_iter_binning = new_region_set
        
        # Update segment data to catch changes in endpoint type
        self.data_manager.update_segments(self.n_iter,self.segments.values())
        
            
    def prepare_new_segments(self):
        self.invoke_callbacks(self.prepare_new_segments)
        self.data_manager.prepare_iteration(self.n_iter+1, self.next_iter_segments)
        
    def run(self):   
        run_starttime = time.time()
        max_walltime = wemd.rc.config.get_interval('limits.max_wallclock', default=None, type=float)
        if max_walltime:
            run_killtime = run_starttime + max_walltime
            wemd.rc.pstatus('Maximum wallclock time: %s' % timedelta(seconds=max_walltime or 0))
        else:
            run_killtime = None
        
        self.n_iter = self.data_manager.current_iteration    
        max_iter = wemd.rc.config.get_int('limits.max_iterations', self.n_iter+1)

        iter_elapsed = 0
        while self.n_iter <= max_iter:
            
            if max_walltime and time.time() + 1.1*iter_elapsed >= run_killtime:
                wemd.rc.pstatus('Iteration {:d} would require more than the allotted time. Ending run.'
                                .format(self.n_iter))
                return
            
            try:
                iter_start_time = time.time()
                
                wemd.rc.pstatus('\n%s' % time.asctime())
                wemd.rc.pstatus('Iteration %d (%d requested)' % (self.n_iter, max_iter))
                                
                self.prepare_iteration()
                wemd.rc.pflush()
                
                self.pre_propagation()
                self.propagate()
                wemd.rc.pflush()
                self.check_propagation()
                wemd.rc.pflush()
                self.post_propagation()
                
                wemd.rc.pflush()
                self.pre_we()
                self.run_we()
                self.post_we()
                wemd.rc.pflush()
                
                self.prepare_new_segments()
                
                self.finalize_iteration()
                
                iter_elapsed = time.time() - iter_start_time
                iter_summary = self.data_manager.get_iter_summary()
                iter_summary['walltime'] += iter_elapsed
                iter_summary['cputime'] = sum(segment.cputime for segment in self.segments.itervalues())
                self.data_manager.update_iter_summary(iter_summary)
    
                self.n_iter += 1
                self.data_manager.current_iteration += 1
                wemd.rc.pflush()
            finally:
                self.data_manager.flush_backing()
                
        wemd.rc.pstatus('\n%s' % time.asctime())
        wemd.rc.pstatus('WEMD run complete.')
        
            
        
    # The functions prepare_run(), finalize_run(), run(), and shutdown() are
    # designed to be called by scripts which are actually performing runs.
    # Specifically, prepare_run() and finalize_run() define the order in which
    # various hooks are called.
    def prepare_run(self):
        '''Prepare a new run.'''
        self.data_manager.prepare_run()
        self.system.prepare_run()
        self.invoke_callbacks(self.prepare_run)
    
    def finalize_run(self):
        '''Perform cleanup at the normal end of a run'''
        self.invoke_callbacks(self.finalize_run)
        self.system.finalize_run()
        self.data_manager.finalize_run()
        
    def pre_propagation(self):
        self.invoke_callbacks(self.pre_propagation)
        
    def post_propagation(self):
        self.invoke_callbacks(self.post_propagation)
        
    def pre_we(self):
        self.invoke_callbacks(self.pre_we)
    
    def post_we(self):
        self.invoke_callbacks(self.post_we)
