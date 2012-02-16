from __future__ import division; __metaclass__ = type

import sys, time, operator, math, numpy, re, random
from itertools import izip
from datetime import timedelta
import logging
log = logging.getLogger(__name__)

import wemd
from wemd.states import BasisState, InitialState, TargetState
from wemd.util import extloader
from wemd import Segment
from wemd.util.miscfn import vgetattr

from wemd.work_managers import ops as wm_ops

EPS = numpy.finfo(numpy.float64).eps

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
                                     self.finalize_iteration, self.prepare_we,
                                     self.prepare_new_segments))
        self._callbacks_by_name = {fn.__name__: fn for fn in self._valid_callbacks}
        self._n_propagated = 0
        
        self.do_gen_istates = wemd.rc.config.get_bool('system.gen_istates', False)
        
        # Per-iteration variables
        self.n_iter = None                  # current iteration
        self.target_states = None           # TargetStates valid at this iteration
        self.target_state_bins = None       # Bins (in self.final_binning) associated with each TargetState
        self.basis_states = None            # BasisStates valid at this iteration
        self.initial_states = None          # InitialStates used in this iteration
        self.spare_initial_states = None    # InitialStates available for use next iteration
        self.assigned_initial_states = None # InitialStates that were available or generated in this iteration but mapped to 
                                            # segments recycled during this iteration 
        
        self.initial_binning = None         # Binning of segments at beginning of this iteration
        self.final_binning = None           # Binning of segments at the end of this iteration
        self.segments = None                # Mapping of seg_id to segment for all segments in this iteration
        self.completed_segments = None      # Mapping of seg_id to segment for all completed segments in this iteration
        self.incomplete_segments = None     # Mapping of seg_id to segment for all incomplete segments in this iteration
        self.to_recycle = None              # Mapping of seg_id to segment for all completed segments to be recycled
        
        self.next_iter_segments = None      # List of segments to be created for the next iteration
        
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
                
    def _invoke_callbacks(self, hook, *args, **kwargs):
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
        
        log.debug('norm is {!r}'.format(norm))
        assert abs(1 - norm) < EPS*len(segments)
        
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


    def initialize_simulation(self, basis_states, target_states, segs_per_state=1, run_we=True):
        '''Initialize a new weighted ensemble simulation, taking ``segs_per_state`` initial
        states from each of the given ``basis_states``.
        
        ``w_init`` is the forward-facing version of this function'''
        
        
        data_manager = self.data_manager
        work_manager = self.work_manager
        propagator = self.propagator
        pstatus = wemd.rc.pstatus
        system = self.system
        
        
        
        pstatus('Calculating progress coordinate values for basis states.')
        futures = [work_manager.submit(wm_ops.get_pcoord, self.propagator, basis_state)
                   for basis_state in basis_states]
        fmap = {future: i for (i, future) in enumerate(futures)}
        for future in work_manager.as_completed(futures): 
            basis_states[fmap[future]].pcoord = future.get_result().pcoord
            
        pstatus('Creating HDF5 file {!r}'.format(self.data_manager.we_h5filename))
        data_manager.prepare_backing()
        data_manager.save_basis_states(basis_states)
        data_manager.save_target_states(target_states)
        
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
        
        pstatus('{:d} target state(s) present'.format(len(target_states)), end='')    
        if wemd.rc.verbose_mode and target_states:
            pstatus(':')
            pstatus('{:6s}    {:12s}    {}'.format('ID', 'Label', 'Progress Coordinate'))
            for target_state in target_states:
                pstatus('{:<6d}    {:12s}    {}'
                                 .format(target_state.state_id, target_state.label, ','.join(map(str,target_state.pcoord))))
        pstatus()

        pstatus('Preparing initial states')
        segments = []
        initial_states = []
        if self.do_gen_istates:
            istate_type = InitialState.ISTATE_TYPE_GENERATED
        else:
            istate_type = InitialState.ISTATE_TYPE_BASIS
            
        for basis_state in basis_states:
            for iseg in xrange(segs_per_state):
                initial_state = data_manager.create_initial_states(0,1)[0]
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
                    
                
        data_manager.update_initial_states(initial_states, n_iter=0)
        
        if run_we:
            # TODO: what to do if something winds up in a recycling region?
            # At the moment, fail with an exception
            region_set = system.new_region_set()
            new_region_set = self.we_driver.run_we(region_set, segments)
        else:
            new_region_set = system.new_region_set()
            new_region_set.assign_to_bins(segments, key=Segment.initial_pcoord)
        
        all_bins = new_region_set.get_all_bins()
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
        data_manager.update_initial_states(initial_states, n_iter=0)
                    
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
        self.initial_binning = self.system.new_region_set()
        self.final_binning = self.system.new_region_set()
        self.target_states = self.data_manager.get_target_states(self.n_iter)
        self.target_state_bins = self.final_binning.map_to_bins([target_state.pcoord for target_state in self.target_states])
        log.debug('target_state_bins={!r}'.format(self.target_state_bins))
        self.basis_states = self.data_manager.get_basis_states(self.n_iter)
        
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
        wemd.rc.pstatus('{:d} segments remain in iteration {:d}'.format(len(incomplete_segments), self.n_iter))
        
        self.initial_binning.assign_to_bins(segments.itervalues(), key=Segment.initial_pcoord)
        self.report_bin_statistics(self.initial_binning, save_summary=True)
        
        self.to_recycle = {}
        if completed_segments:
            self.final_binning.assign_to_bins(completed_segments.values(), key=Segment.final_pcoord)
            
            for bin, target_state in izip(self.target_state_bins, self.target_states):
                log.debug('bin={!r}, target_state={!r}'.format(bin,target_state))
                if len(bin):
                    log.debug('{:d} replicas in target state {!r}'.format(len(bin), target_state))
                    self.to_recycle.update({segment.seg_id: segment for segment in bin})
        log.debug('{:d} replicas in target states'.format(len(self.to_recycle)))
        
        self.initial_states = self.data_manager.get_segment_initial_states(segments.values())
        
        self.assigned_initial_states = set()
        self.spare_initial_states = set(self.data_manager.get_unused_initial_states())
        log.debug('{:d} unused initial states found'.format(len(self.spare_initial_states)))
        
        log.debug('dispatching propagator prep_iter to work manager')
        self.work_manager.submit(wm_ops.prep_iter, self.propagator, self.n_iter, segments).get_result()
        
    def finalize_iteration(self, n_iter, segments):
        '''Perform customized processing/cleanup on just-completed segments at the end of an iteration'''
        log.debug('finalizing iteration {:d}'.format(self.n_iter))
        self.work_manager.submit(wm_ops.post_iter, self.propagator, self.n_iter, self.segments.values()).get_result()


    def ensure_istates(self):
        '''Ensure that as many unused initial states are available as there are currently 
        particles to be recycled.  This will select a basis state and (depending on whether
        initial state generation is on) will either immediately make an initial state available
        or dispatch a task to the work manager to generate a state from the basis state.
        
        Returns a set of futures corresponding to initial state generation tasks (if any)
        '''
        
        # Move as many spare states as necessary from the spare pool to the assigned pool
        while len(self.to_recycle) > len(self.assigned_initial_states) and self.spare_initial_states:
            self.assigned_initial_states.add(self.spare_initial_states.pop())
            
        # If we still need more, we need to generate them
        n_needed = len(self.to_recycle) - len(self.assigned_initial_states)
        if n_needed < 1:
            return set()

        if self.do_gen_istates:
            istate_type = InitialState.ISTATE_TYPE_GENERATED
        else:
            istate_type = InitialState.ISTATE_TYPE_BASIS
            
        futures = set()
        cumul_probs = numpy.add.accumulate([bstate.probability for bstate in self.basis_states])
        for inew in xrange(n_needed):
            # Select a basis state according to its weight
            ibstate = numpy.digitize([random.random()], cumul_probs)
            basis_state = self.basis_states[ibstate]
            initial_state = self.data_manager.create_initial_states(self.n_iter,1)[0]
            initial_state.basis_state_id = basis_state.state_id
            initial_state.istate_type = istate_type
            initial_state.basis_state = basis_state
            
            if self.do_gen_istates:
                # We need to dispatch a task to the work manager to generate a new state
                futures.add(self.work_manager.submit(wm_ops.gen_istate, self.propagator, initial_state.basis_state, initial_state))
            else:
                # We can immediately generate a new state, mark it as ready, and save it
                initial_state.pcoord = basis_state.pcoord
                initial_state.istate_status = InitialState.ISTATE_STATUS_PREPARED
                self.data_manager.update_initial_states([initial_state])
                self.assigned_initial_states.add(initial_state)
                self.initial_states[initial_state.state_id: initial_state]
                
        assert len(self.to_recycle) == len(self.assigned_initial_states) + len(futures)
        return futures
                        
    def propagate(self):
        segments = self.incomplete_segments.values()
        log.debug('iteration {:d}: propagating {:d} segments'.format(self.n_iter, len(segments)))
        futures = set()
        segment_futures = set()
        istate_gen_futures = set()
                
        for segment in segments:
            pbstates, pistates = wemd.states.pare_basis_initial_states(self.basis_states, self.initial_states, [segment])
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
                self.data_manager.update_segments(self.n_iter, incoming)
                self.data_manager.flush_backing()
                
                self.segments.update({segment.seg_id: segment for segment in incoming})
                self.completed_segments.update({segment.seg_id: segment for segment in incoming})
                
                for segment in incoming:
                    final_bin = self.final_binning.map_to_bins([segment.pcoord[-1]])[0]
                    log.debug('incoming segment {!r} mapped to bin {!r}'.format(segment, final_bin))
                    final_bin.add(segment)
                    if final_bin in self.target_state_bins:
                        # This particle must be recycled
                        self.to_recycle[segment.seg_id: segment]
                        new_futures = self.ensure_istates()
                        istate_gen_futures.update(new_futures)
                        futures.update(new_futures)
                
            elif future in istate_gen_futures:
                istate_gen_futures.remove(future)
                _basis_state, initial_state = future.get_result()
                initial_state.istate_status = InitialState.ISTATE_STATUS_PREPARED
                self.data_manager.update_initial_states([initial_state])
                self.assigned_initial_states.add(initial_state)
                self.initial_states[initial_state.state_id: initial_state]
                self.data_manager.flush_backing()
            else:
                raise AssertionError('untracked future')                    
                    
        log.debug('done with propagation')

    def check_propagation(self):
        failed_segments = [segment for segment in self.segments.itervalues() if segment.status != Segment.SEG_STATUS_COMPLETE]
        
        if failed_segments:
            failed_ids = '  \n'.join(str(segment.seg_id) for segment in failed_segments)
            log.error('propagation failed for {:d} segment(s):\n{}'.format(len(failed_segments), failed_ids))
            raise PropagationError('propagation failed for {:d} segments'.format(len(failed_segments)))
        else:
            log.debug('propagation complete for iteration {:d}'.format(self.n_iter))
            
        failed_istates = [istate for istate in self.assigned_initial_states 
                          if istate.istate_status != InitialState.ISTATE_STATUS_PREPARED]
        if failed_istates:
            failed_ids = '  \n'.join(str(istate.state_id) for istate in failed_istates)
            log.error('initial state generation failed for {:d} states:\n{}'.format(len(failed_istates, failed_ids)))
            raise PropagationError('initial state generation failed for {:d} states'.format(len(failed_istates)))
        else:
            log.debug('initial state generation complete for iteration {:d}'.format(self.n_iter))


    def run_we(self):
        '''Run the weighted ensemble algorithm based on the binning in self.final_bins and
        the recycled particles in self.to_recycle.'''


        # assign initial states to particles being recycled
        init_segs = []
        used_initial_states = set()
        for segment, initial_state in izip(self.to_recycle.itervalues(), self.assigned_initial_states):
            new_segment = Segment(n_iter=self.n_iter+1, parent_id=-(initial_state.state_id+1),
                                  weight=segment.weight, pcoord=self.system.new_pcoord_array(),
                                  wtg_parent_ids=[-(initial_state.state_id+1)],)
            new_segment.pcoord[0] = initial_state.pcoord
            init_segs.append(new_segment)
                        
        new_region_set = self.we_driver.run_we(self.final_binning, init_segs)
        new_segments = list(new_region_set.particles)
        
        for segment in new_segments:
            segment.n_iter = self.n_iter+1
            segment.status = Segment.SEG_STATUS_PREPARED
            if segment.initpoint_type == Segment.SEG_INITPOINT_NEWTRAJ:
                initial_state = self.initial_states[segment.initial_state_id]
                initial_state.iter_used = self.n_iter+1
                used_initial_states.add(initial_state)            
                        
        self.data_manager.prepare_iteration(self.n_iter+1, new_segments)
        if used_initial_states:
            self.data_manager.update_initial_states(used_initial_states)
        
                
                
            
        

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

        
        while self.n_iter < max_iter:
            iter_start_time = time.time()
            iter_summary = self.data_manager.get_iter_summary()
            self.prepare_iteration()
            self.propagate()
            self.check_propagation()
            self.run_we()
            
            
            iter_elapsed = time.time() - iter_start_time
            iter_summary['walltime'] += iter_elapsed
            iter_summary['cputime'] = sum(segment.cputime for segment in self.segments.itervalues())
            self.data_manager.update_iter_summary(iter_summary)

            self.n_iter += 1
            self.data_manager.current_iteration += 1        
            
    def _propagate(self, n_iter, segments):            
        wemd.rc.pflush()
        log.debug('propagating {:d} segments'.format(len(segments)))
        futures = []
        for segment in segments:
            future = self.work_manager.submit(wm_ops.propagate, self.propagator, [segment])
            log.debug('segment {:d} is task {!r}'.format(segment.seg_id, future.task_id))
            futures.append(future)
        completed = []
        for future in self.work_manager.as_completed(futures):
            incoming = future.get_result()
            self._n_propagated += 1
            log.debug('recording results for {!r}, {:d} for this run'.format(incoming, self._n_propagated))
            self.data_manager.update_segments(n_iter, incoming)
            self.data_manager.flush_backing()
            completed.extend(incoming)
        log.debug('done with propagation')
        assert len(completed) == len(segments)
        return completed        
        
    # The functions prepare_run(), finalize_run(), run(), and shutdown() are
    # designed to be called by scripts which are actually performing runs.
    # Specifically, prepare_run() and finalize_run() define the order in which
    # various hooks are called.
    def prepare_run(self):
        '''Prepare a new run.'''
        self.data_manager.prepare_run()
        self.system.prepare_run()
        self._invoke_callbacks(self.prepare_run)
    
    def finalize_run(self):
        '''Perform cleanup at the normal end of a run'''
        self._invoke_callbacks(self.finalize_run)
        self.system.finalize_run()
        self.data_manager.finalize_run()
    
    """    
    def prepare_iteration(self, n_iter, segments, partial):
        '''Perform customized processing/setup prior to propagation. Argument ``partial`` (True or False)
        indicates whether this is a partially-complete iteration (i.e. a restart)'''
        self.work_manager.submit(wm_ops.prep_iter, self.propagator, n_iter, segments).get_result()
        
    def prepare_propagation(self, n_iter, segments):
        '''Prepare to propagate a group of segments'''
        pass
    
    def finalize_propagation(self, n_iter, segments):
        '''Clean up after propagation'''
        pass
    
    def finalize_iteration(self, n_iter, segments):
        '''Perform customized processing/cleanup on just-completed segments at the end of an iteration'''
        self.work_manager.submit(wm_ops.post_iter, self.propagator, n_iter, segments).get_result()
    """    
    def prepare_we(self, n_iter, segments):
        '''Perform customized processing after propagation and before WE'''
        self._invoke_callbacks(self.prepare_we, n_iter, segments)
    
    def prepare_new_segments(self, n_iter, segments, next_iter_segments):
        self._invoke_callbacks(self.prepare_new_segments, n_iter, segments, next_iter_segments)
    
    def shutdown(self, exit_code = 0):
        '''Shut down a simulation in an orderly way'''
        self.work_manager.shutdown(exit_code)
        
                            
    def old_run(self):
        """Begin (or continue) running a simulation.  Must only be called in master processes.
        """
                
        # Set up internal timing
        run_starttime = time.time()
        max_walltime = wemd.rc.config.get_interval('limits.max_wallclock', default=None, type=float)
        if max_walltime:
            run_killtime = run_starttime + max_walltime
            wemd.rc.pstatus('Maximum wallclock time: %s' % timedelta(seconds=max_walltime or 0))
        else:
            run_killtime = None
        iteration_elapsed = 0
        
        # Wrap everything in try/finally to ensure that the data manager gets things written
        # correctly (should work without this, but the h5py documentation is a little hazy
        # on how wise it is)
        try:            
            # Get segments
            n_iter = self.data_manager.current_iteration
            max_iter = wemd.rc.config.get_int('limits.max_iterations', n_iter+1)
                         
            # Guaranteed ordering by seg_id, so segments[seg_id] works for any valid seg_id for this iteration
            segments = self.data_manager.get_segments(n_iter)
            while n_iter <= max_iter:
                log.debug('beginning iteration {:d}'.format(n_iter))
                
                # Check to see if we will exceed the allowable wallclock time
                if max_walltime and time.time() + 1.1*iteration_elapsed >= run_killtime:
                    wemd.rc.pstatus('Iteration %d would require more than the allotted time. Ending run.'
                                             % n_iter)
                    self.shutdown(0)
                    return
                
                # Record/report time and iteration number    
                iteration_starttime = time.time()
                wemd.rc.pstatus('\n%s' % time.asctime())
                wemd.rc.pstatus('Iteration %d (%d requested)' % (n_iter, max_iter))
                
                # Check probabilities and report
                seg_weights = vgetattr('weight', segments, numpy.float64)
                assert not (seg_weights == 0).any()
                norm = seg_weights.sum()
                wemd.rc.pstatus('norm: %g, error in norm: %g' % (norm, norm-1))
                
                # Determine how many segments remain in this iteration
                segs_to_run = [segment for segment in segments if segment.status != Segment.SEG_STATUS_COMPLETE]
                wemd.rc.pstatus('%d of %d segments remain in iteration %d' % (len(segs_to_run), len(segments), n_iter))
                partial_iteration = (len(segs_to_run) != len(segments))
                
                self.prepare_iteration(n_iter, segments, partial_iteration)
                
                if not partial_iteration:
                    # First run within this iteration; store bin distribution statistics in HDF5
                                
                    # All bins will be empty if this is the first iteration in this run (i.e. this [Unix] process)
                    # If so, perform binning to get statistics to store
                    bins = self.system.region_set.get_all_bins()
                    bin_counts = vgetattr('count', bins, numpy.uint)
                    target_counts = vgetattr('target_count', bins, numpy.uint)
                    if (bin_counts == 0).all():                
                        log.info('initial iteration for this run; binning on segment initial points')
                        initial_points = [segment.pcoord[0] for segment in segments]
                        for (segment, bin) in izip(segments, self.system.region_set.map_to_bins(initial_points)):
                            bin.add(segment)
                        bin_counts = vgetattr('count', bins, numpy.uint)
    
                    # Do not include bins with target count zero (e.g. sinks, never-filled bins) in the (non)empty bins statistics
                    n_active_bins = len(target_counts[target_counts!=0])
                    seg_probs  = vgetattr('weight', segments, numpy.float64)
                    bin_probs  = vgetattr('weight', bins, numpy.float64)
                    
                    assert abs(1 - seg_probs.sum()) < EPS*len(segments)
                    
                    min_seg_prob = seg_probs[seg_probs!=0].min()
                    max_seg_prob = seg_probs.max()
                    seg_drange   = math.log(max_seg_prob/min_seg_prob)
                    min_bin_prob = bin_probs[bin_probs!=0].min()
                    max_bin_prob = bin_probs.max()
                    bin_drange = math.log(max_bin_prob/min_bin_prob)
                    n_pop = len(bin_counts[bin_counts!=0])
                    wemd.rc.pstatus('{:d} of {:d} ({:%}) active bins are populated'.format(n_pop, n_active_bins, 
                                                                                                      n_pop/n_active_bins))
                    wemd.rc.pstatus('per-bin minimum non-zero probability:       {:g}'.format(min_bin_prob))
                    wemd.rc.pstatus('per-bin maximum probability:                {:g}'.format(max_bin_prob))
                    wemd.rc.pstatus('per-bin probability dynamic range (kT):     {:g}'.format(bin_drange))
                    wemd.rc.pstatus('per-segment minimum non-zero probability:   {:g}'.format(min_seg_prob))
                    wemd.rc.pstatus('per-segment maximum non-zero probability:   {:g}'.format(max_seg_prob))
                    wemd.rc.pstatus('per-segment probability dynamic range (kT): {:g}'.format(seg_drange))
                                        
                    iter_summary = self.data_manager.get_iter_summary(n_iter)
                    iter_summary['n_particles'] = len(segments)
                    iter_summary['norm'] = norm
                    iter_summary['min_bin_prob'] = min_bin_prob
                    iter_summary['max_bin_prob'] = max_bin_prob
                    iter_summary['min_seg_prob'] = min_seg_prob
                    iter_summary['max_seg_prob'] = max_seg_prob
                    self.data_manager.update_iter_summary(n_iter, iter_summary)
                
                # Allow the user to run only one segment to aid in debugging propagators
                log.debug('propagating iteration {:d}'.format(n_iter))
                if wemd.rc.config.get('args.only_one_segment',False):
                    log.info('propagating only one segment')
                    segments_run = self.propagate(n_iter, segs_to_run[0:1])                
                    if len(segs_to_run) > 1:
                        # Exit cleanly after one segment, but only if we haven't finished the last segment in an iteration
                        # In the latter case, go ahead and run WE
                        break
                else:
                    segments_run = self.propagate(n_iter, segs_to_run)
                del segs_to_run
                log.debug('propagation complete')
                
                # Collect segments into one place again
                segs_by_id = {segment.seg_id: segment for segment in segments}
                for segment in segments_run:
                    segs_by_id[segment.seg_id] = segment
                segments = sorted(segs_by_id.itervalues(), key=operator.attrgetter('seg_id'))
                assert all(segment.seg_id == i for i, segment in enumerate(segments))
                                
                # Check to ensure that all segments have been propagated                    
                failed_segments = [segment for segment in segments if segment.status != Segment.SEG_STATUS_COMPLETE]
                if failed_segments:
                    wemd.rc.pstatus('Propagation FAILED for %d segments:' % len(failed_segments))
                    for failed_segment in failed_segments:
                        wemd.rc.pstatus('  %d' % failed_segment.seg_id)
                    raise RuntimeError('propagation failed for %d segments' % len(failed_segments))
                
                region_set = self.system.region_set
                region_set.clear()
                bins = region_set.get_all_bins()
                n_bins = len(bins)
                n_segs = len(segments)
                pcoord_len = self.system.pcoord_len
                pcoord_dtype = self.system.pcoord_dtype
                pcoord_ndim = self.system.pcoord_ndim
                
                # Progress coordinates for all segments (seg_id, time in segment, dimension) -> progress coordinate vector
                pcoords = numpy.empty((n_segs, pcoord_len, pcoord_ndim), dtype=pcoord_dtype)
                # Assignments to bins at each time point (seg_id, time in segment) -> bin index
                assignments = numpy.empty((n_segs, pcoord_len), numpy.uint32)
                # Instantaneous bin populations (time point, bin index) -> population
                populations = numpy.zeros((pcoord_len, n_bins), numpy.float64)
                # n_trans[i,j] = number of transitions from bin i at the beginning of this iteration to bin j at the end
                n_trans = numpy.zeros((n_bins,n_bins), numpy.uint32)
                # fluxes[i,j] = probability flow from bin i at the beginning of this iteration to bin j at the end
                fluxes = numpy.zeros((n_bins,n_bins), numpy.float64)
                # rate[i,j] = measured rate from bin i to bin j
                rates = numpy.zeros((n_bins,n_bins), numpy.float64)

                # Pull progress coordinate data from this iteration into one big array
                for (seg_id, segment) in enumerate(segments):
                    pcoords[seg_id,...] = segment.pcoord[...]
                    
                # Assign progress coordinates to bins
                for i in xrange(0, pcoord_len):
                    assignments[:,i] = region_set.map_to_all_indices(pcoords[:,i])
                    
                # Calculate instantaneous bin probabilities
                for (seg_id, segment) in enumerate(segments):
                    weight = segment.weight
                    for i in xrange(0, pcoord_len):
                        populations[i, assignments[seg_id,i]] += weight
                    
                # Calculate number of transitions, fluxes, and rates
                # Simultaneously put segments into their respective bins
                for (segment, init_assignment, final_assignment) in izip(segments, assignments[:,0], assignments[:,-1]):
                    n_trans[init_assignment, final_assignment] += 1
                    
                    # Flux is in units of tau^-1
                    fluxes[init_assignment, final_assignment] += segment.weight
                    
                    # rate[i,j] = flux[i,j] / population[i] (at beginning of iteration)
                    # note the implicit loop (broadcast) over j
                    for i in xrange(0, n_bins):
                        rates[i,:] = fluxes[i,:] / populations[0,i] if populations[0,i] > 0 else 0
                        
                    bins[final_assignment].add(segment)
                        
                self.data_manager.write_bin_data(n_iter, assignments, populations, n_trans, fluxes, rates)


                
                seg_probs  = vgetattr('weight', segments, numpy.float64)
                bin_probs  = vgetattr('weight', bins, numpy.float64)
                assert abs(seg_probs.sum() - 1) < EPS*len(segments)
                                        
                self.prepare_we(n_iter, segments)

                seg_probs  = vgetattr('weight', segments, numpy.float64)
                bin_probs  = vgetattr('weight', bins, numpy.float64)
                log.debug('norm of seg probs prior to WE: {:20.16f}'.format(seg_probs.sum()))
                log.debug('norm of bin probs prior to WE: {:20.16f}'.format(bin_probs.sum()))
                assert abs(seg_probs.sum() - 1) < EPS*len(segments)
                
                # Run the weighted ensemble algorithm
                next_iter_segments = self.we_driver.run_we(segments, self.system.region_set)
                
                next_seg_probs  = vgetattr('weight', next_iter_segments, numpy.float64)
                bin_probs = vgetattr('weight', bins, numpy.float64)
                log.debug('norm of seg probs after WE: {:20.16f}'.format(next_seg_probs.sum()))
                log.debug('norm of bin probs after WE: {:20.16f}'.format(bin_probs.sum()))
                
                assert abs(next_seg_probs.sum() - 1) < EPS*len(next_iter_segments)
                                
                # ``segments`` is still set of segments in just-completed iteration, 
                # but now ``region_set`` contains next iteration's segments
                                
                # Report recycling statistics
                n_recycled = sum(dest.count for dest in self.we_driver.recycle_to)
                P_recycled = sum(dest.weight for dest in self.we_driver.recycle_to)                            
                if n_recycled > 0:
                    P_recycled = sum(dest.weight for dest in self.we_driver.recycle_to)
                    wemd.rc.pstatus('%d particles (%g probability) recycled' % (n_recycled, P_recycled))
                    log.debug('from: %r' % self.we_driver.recycle_from)
                    log.debug('to:   %r' % self.we_driver.recycle_to)
                    for isrc,src in enumerate(self.we_driver.recycle_from):
                        wemd.rc.pstatus("  {0.count:d} ({0.weight:g} probability) from region {1:d} '{2}'"\
                                                 .format(src,isrc,self.system.target_states[isrc].label))
                        
                    for idest, dest in enumerate(self.we_driver.recycle_to):
                        wemd.rc.pstatus(("  {0.count:d} ({0.weight:g} probability,"
                                                  +" {1:%} of recycled particles) to state '{2}'")\
                                                 .format(dest, dest.count/n_recycled,
                                                         self.system.initial_states[idest].label))
                else:
                    wemd.rc.pstatus('0 particles recycled')
                    
                wemd.rc.pstatus('%d trajectory segments in next iteration' % len(next_iter_segments))
                
                # Store recycling information in HDF5, but only if we actually have targets
                iter_summary = self.data_manager.get_iter_summary(n_iter)                
                if self.system.target_states:    
                    iter_summary['target_flux'] = P_recycled
                    iter_summary['target_hits'] = n_recycled
                    self.data_manager.write_recycling_data(n_iter, self.we_driver.recycle_from)
                
                # Update current segments; their endpoint types were modified by WE                                
                self.data_manager.update_segments(n_iter, segments)                
                self.data_manager.flush_backing()
                                
                log.debug('work manager finalizing iteration')
                self.finalize_iteration(n_iter, segments)
                # For use in debugging
                seg_debug_fmt = '\n            '.join(['created segment {segment!r}', 
                                                 'n_iter: {segment.n_iter}, seg_id: {segment.seg_id}, weight: {segment.weight}',
                                                 'n_parents: {segment.n_parents}, p_parent_id: {segment.p_parent_id}, '
                                                 +'parent_ids: {segment.parent_ids!r}',
                                                 ' pcoord[0]: {init_pcoord!r}',
                                                 'pcoord[-1]: {final_pcoord!r}'
                                                 ])
            
                if log.isEnabledFor(logging.DEBUG):
                    log.debug('new segments follow')
                for new_segment in next_iter_segments:
                    assert new_segment.weight is not None
                    assert new_segment.p_parent_id is not None
                    assert new_segment.parent_ids is not None and new_segment.p_parent_id in new_segment.parent_ids
                    new_segment.n_iter = n_iter+1
                    new_segment.status = Segment.SEG_STATUS_PREPARED
                    new_segment.n_parents = len(new_segment.parent_ids)
                
                    if log.isEnabledFor(logging.DEBUG):
                        log.debug(seg_debug_fmt.format(segment = new_segment, init_pcoord = segment.pcoord[0],
                                                       final_pcoord = segment.pcoord[-1]))
                
                # Do any modification necessary prior to committing the next iteration's segments to storage
                self.prepare_new_segments(n_iter, segments, next_iter_segments)
                        
                # Create new iteration group in HDF5            
                self.data_manager.prepare_iteration(n_iter+1, next_iter_segments)
                
                # Store timing information 
                iter_walltime = time.time() - iteration_starttime
                iter_cputime = sum(segment.cputime or 0.0 for segment in segments)
                
                iter_summary['walltime'] += iter_walltime
                iter_summary['cputime'] = iter_cputime
                self.data_manager.update_iter_summary(n_iter, iter_summary)
                
                try:
                    #This may give NaN if starting a truncated simulation
                    walltime = timedelta(seconds=float(iter_summary['walltime']))
                except ValueError:
                    walltime = 0.0 
                
                try:
                    cputime = timedelta(seconds=float(iter_summary['cputime']))
                except ValueError:
                    cputime = 0.0                     
                    
                wemd.rc.pstatus('Iteration wallclock: {0!s}, cputime: {1!s}'\
                                          .format(walltime,
                                                  cputime))
                
                # Update HDF5 and flush the status output buffer in preparation for
                # the next iteration
                n_iter += 1
                self.data_manager.current_iteration = n_iter
                self.data_manager.flush_backing()
                wemd.rc.pflush()
                
                # Update segments list for next iteration
                segments = next_iter_segments            
                
            # end propagation/WE loop
        finally:
            # Regardless of what happens (end of loop/exit/exception), make sure run data
            # is saved
            self.data_manager.flush_backing()
            self.data_manager.close_backing()
        
        wemd.rc.pstatus('\n%s' % time.asctime())
        wemd.rc.pstatus('WEMD run complete.')
    # end WESimManager.run()
# end SimManager
    
