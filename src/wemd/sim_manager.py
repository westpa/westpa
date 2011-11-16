from __future__ import division; __metaclass__ = type

import sys, time, operator, math, numpy, re
from itertools import izip, imap
from datetime import timedelta
import logging
log = logging.getLogger(__name__)

import wemd
from wemd.util import extloader
from wemd.util.rtracker import ResourceTracker, ResourceUsage
from wemd import Segment
from wemd.util.miscfn import vgetattr

def wm_prep_iter(propagator, n_iter):
    propagator.pre_iter(n_iter)
    
def wm_post_iter(propagator, n_iter):
    propagator.post_iter(n_iter)
    
def wm_propagate(propagator, segments):
    '''Propagate the given segments with the given propagator. This has to be a top-level
    function for the current incarnation of the work manager.'''
    outgoing_ids = [segment.seg_id for segment in segments]
    incoming_segments = {segment.seg_id: segment for segment in propagator.propagate(segments)}
    log.debug('received {!r}'.format(incoming_segments))
    return [incoming_segments[seg_id] for seg_id in outgoing_ids]


class WESimManager:
    def __init__(self):        
        self.data_manager = None
        self.we_driver = None
        self.work_manager = None        
        self.propagator = None
        
        self.system = None
                        
        self.status_stream = sys.stdout
        
        self.rtracker = ResourceTracker()
                
        # A table of function -> list of (priority, name, callback) tuples
        self._callback_table = {}
        self._valid_callbacks = set((self.prepare_run, self.finalize_run,
                                     self.finalize_iteration, self.prepare_we,
                                     self.prepare_new_segments))
        self._callbacks_by_name = {fn.__name__: fn for fn in self._valid_callbacks}
        
    def register_callback(self, hook, function, priority=0):
        '''Registers a callback to execute during the given ``hook`` into the simulation loop. The optional
        priority is used to order when the function is called relative to other registered callbacks.'''
        
        if hook not in self._valid_callbacks:
            try:
                hook = self._callbacks_by_name[hook]
            except KeyError:
                raise KeyError('invalid hook {!r}'.format(hook))
            
        try:
            self._callback_table[hook].add([(priority,function.__name__,function)])
        except KeyError:
            self._callback_table[hook] = set([(priority,function.__name__,function)])
        
        log.debug('registered callback {!r} for hook {!r}'.format(function, hook))
        
    def _invoke_callbacks(self, hook, *args, **kwargs):
        callbacks = self._callback_table.get(hook, [])
        sorted_callbacks = list(sorted(callbacks))
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
        
    def flush_status(self):
        try:
            self.status_stream.flush()
        except AttributeError:
            pass
    
    def propagate(self, n_iter, segments):
        # Propagate segments            
        self.rtracker.begin('propagation')
        self.flush_status()

        futures = [self.work_manager.submit(wm_propagate, self.propagator, [segment]) for segment in segments]
        completed = []
        #self.work_manager.wait_all(futures)
        for future in self.work_manager.as_completed(futures):
        #for future in futures:
            incoming = future.get_result()
            self.data_manager.update_segments(n_iter, incoming)
            completed.extend(incoming)
        
        self.rtracker.end('propagation')
        return completed
        
        

    # The functions prepare_run(), finalize_run(), run(), and shutdown() are
    # designed to be called by scripts which are actually performing runs.
    # Specifically, prepare_run() and finalize_run() define the order in which
    # various hooks are called.
    def prepare_run(self):
        '''Prepare a new run.'''
        self.propagator.system = self.system
        self.work_manager.propagator = self.propagator
        self.data_manager.system = self.system
        self.we_driver.system = self.system
        
        #self.work_manager.prepare_run()
        self.data_manager.prepare_run()
        self.system.prepare_run()
        self._invoke_callbacks(self.prepare_run)
    
    def finalize_run(self):
        '''Perform cleanup at the normal end of a run'''
        self.system.finalize_run()
        self.data_manager.finalize_run()
        #self.work_manager.finalize_run()
        self._invoke_callbacks(self.finalize_run)
        
    def prepare_iteration(self, n_iter, segments, partial):
        '''Perform customized processing/setup prior to propagation. Argument ``partial`` (True or False)
        indicates whether this is a partially-complete iteration (i.e. a restart)'''
        self.work_manager.submit(wm_prep_iter, self.propagator, n_iter).wait()
        
    def prepare_propagation(self, n_iter, segments):
        '''Prepare to propagate a group of segments'''
        pass
    
    def finalize_propagation(self, n_iter, segments):
        '''Clean up after propagation'''
        pass
    
    def finalize_iteration(self, n_iter, segments):
        '''Perform customized processing/cleanup on just-completed segments at the end of an iteration'''
        self.work_manager.submit(wm_post_iter, self.propagator, n_iter).wait()
        
    def prepare_we(self, n_iter, segments):
        '''Perform customized processing after propagation and before WE'''
        self._invoke_callbacks(self.prepare_we, n_iter, segments)
    
    def prepare_new_segments(self, n_iter, segments, next_iter_segments):
        self._invoke_callbacks(self.prepare_new_segments, n_iter, segments, next_iter_segments)
    
    def shutdown(self, exit_code = 0):
        '''Shut down a simulation in an orderly way'''
        self.work_manager.shutdown(exit_code)
                            
    def run(self):
        """Begin (or continue) running a simulation.  Must only be called in master processes.
        """
        # Set up internal timing
        self.rtracker.begin('run')
        run_starttime = time.time()
        max_walltime = wemd.rc.config.get_interval('limits.max_wallclock', default=None, type=float)
        if max_walltime:
            run_killtime = run_starttime + max_walltime
            self.status_stream.write('Maximum wallclock time: %s\n' % timedelta(seconds=max_walltime or 0))
        else:
            run_killtime = None
        iteration_elapsed = 0
        
        # Wrap everything in try/finally to ensure that the data manager gets things written
        # correctly (should work without this, but the h5py documentation is a little hazy
        # on how wise it is)
        try:
            self.data_manager.open_backing()
                        
            # Get segments
            n_iter = self.data_manager.current_iteration
            max_iter = wemd.rc.config.get_int('limits.max_iterations', n_iter+1)
                         
            # Guaranteed ordering by seg_id, so segments[seg_id] works for any valid seg_id for this iteration
            segments = self.data_manager.get_segments(n_iter)
            while n_iter <= max_iter:
                self.rtracker.begin('iteration')
                
                # Check to see if we will exceed the allowable wallclock time
                if max_walltime and time.time() + 1.1*iteration_elapsed >= run_killtime:
                    self.status_stream.write('Iteration %d would require more than the allotted time. Ending run.\n'
                                             % n_iter)
                    self.shutdown(0)
                    return
                
                # Record/report time and iteration number    
                iteration_starttime = time.time()
                self.status_stream.write('\n%s\n' % time.asctime())
                self.status_stream.write('Iteration %d (%d requested)\n' % (n_iter, max_iter))
                
                # Check probabilities and report
                seg_weights = vgetattr('weight', segments, numpy.float64)
                assert not (seg_weights == 0).any()
                norm = seg_weights.sum()
                self.status_stream.write('norm: %g, error in norm: %g\n' % (norm, norm-1))
                
                # Determine how many segments remain in this iteration
                segs_to_run = [segment for segment in segments if segment.status != Segment.SEG_STATUS_COMPLETE]
                self.status_stream.write('%d of %d segments remain in iteration %d\n' % (len(segs_to_run), len(segments), n_iter))
                partial_iteration = bool(not (len(segs_to_run) == len(segments)))
                
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
                    
                    assert abs(1 - seg_probs.sum()) < 1.0e-15*len(segments)
                    
                    min_seg_prob = seg_probs[seg_probs!=0].min()
                    max_seg_prob = seg_probs.max()
                    seg_drange   = math.log(max_seg_prob/min_seg_prob)
                    min_bin_prob = bin_probs[bin_probs!=0].min()
                    max_bin_prob = bin_probs.max()
                    bin_drange = math.log(max_bin_prob/min_bin_prob)
                    n_pop = len(bin_counts[bin_counts!=0])
                    self.status_stream.write('{:d} of {:d} ({:%}) active bins are populated\n'.format(n_pop, n_active_bins, 
                                                                                                      n_pop/n_active_bins))
                    self.status_stream.write('per-bin minimum non-zero probability:       {:g}\n'.format(min_bin_prob))
                    self.status_stream.write('per-bin maximum probability:                {:g}\n'.format(max_bin_prob))
                    self.status_stream.write('per-bin probability dynamic range (kT):     {:g}\n'.format(bin_drange))
                    self.status_stream.write('per-segment minimum non-zero probability:   {:g}\n'.format(min_seg_prob))
                    self.status_stream.write('per-segment maximum non-zero probability:   {:g}\n'.format(max_seg_prob))
                    self.status_stream.write('per-segment probability dynamic range (kT): {:g}\n'.format(seg_drange))
                                        
                    iter_summary = self.data_manager.get_iter_summary(n_iter)
                    iter_summary['n_particles'] = len(segments)
                    iter_summary['norm'] = norm
                    iter_summary['min_bin_prob'] = min_bin_prob
                    iter_summary['max_bin_prob'] = max_bin_prob
                    iter_summary['bin_dyn_range'] = bin_drange
                    iter_summary['min_seg_prob'] = min_seg_prob
                    iter_summary['max_seg_prob'] = max_seg_prob
                    iter_summary['seg_dyn_range'] = seg_drange
                    self.data_manager.update_iter_summary(n_iter, iter_summary)
                    
                
                # Allow the user to run only one segment to aid in debugging propagators
                log.debug('propagating iteration {:d}'.format(n_iter))
                if wemd.rc.config.get('args.only_one_segment',False):
                    log.info('propagating only one segment')
                    segments = self.propagate(n_iter, segs_to_run[0:1])                
                    if len(segs_to_run) > 1:
                        # Exit cleanly after one segment, but only if we haven't finished the last segment in an iteration
                        # In that case, go ahead and run WE
                        break
                else:
                    segments = self.propagate(n_iter, segs_to_run)
                log.debug('propagation complete')
                                
                # Check to ensure that all segments have been propagated                    
                failed_segments = [segment for segment in segments if segment.status != Segment.SEG_STATUS_COMPLETE]
                if failed_segments:
                    self.status_stream.write('Propagation FAILED for %d segments:\n' % len(failed_segments))
                    for failed_segment in failed_segments:
                        self.status_stream.write('  %d\n' % failed_segment.seg_id)
                    raise RuntimeError('propagation failed for %d segments' % len(failed_segments))
                
                # Sort segments by id
                segments = sorted(segments, key=operator.attrgetter('seg_id'))
                assert all(segment.seg_id == i for i, segment in enumerate(segments))
                
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
                assert abs(seg_probs.sum() - 1) < 1.0e-15*len(segments)
                                        
                self.rtracker.begin('we_prep')
                self.prepare_we(n_iter, segments)
                self.rtracker.end('we_prep')

                seg_probs  = vgetattr('weight', segments, numpy.float64)
                bin_probs  = vgetattr('weight', bins, numpy.float64)
                log.debug('norm of seg probs prior to WE: {:20.16f}'.format(seg_probs.sum()))
                log.debug('norm of bin probs prior to WE: {:20.16f}'.format(bin_probs.sum()))
                assert abs(seg_probs.sum() - 1) < 1.0e-15*len(segments)
                
                # Run the weighted ensemble algorithm
                self.rtracker.begin('we_core')
                next_iter_segments = self.we_driver.run_we(segments, self.system.region_set)
                self.rtracker.end('we_core')
                
                next_seg_probs  = vgetattr('weight', next_iter_segments, numpy.float64)
                bin_probs = vgetattr('weight', bins, numpy.float64)
                log.debug('norm of seg probs after WE: {:20.16f}'.format(next_seg_probs.sum()))
                log.debug('norm of bin probs after WE: {:20.16f}'.format(bin_probs.sum()))
                
                assert abs(next_seg_probs.sum() - 1) < 1.0e-15*len(next_iter_segments)
                                
                # ``segments`` is still set of segments in just-completed iteration, 
                # but now ``region_set`` contains next iteration's segments
    
                self.rtracker.begin('we_postprocess')
                                
                # Report recycling statistics
                n_recycled = sum(dest.count for dest in self.we_driver.recycle_to)
                P_recycled = sum(dest.weight for dest in self.we_driver.recycle_to)                            
                if n_recycled > 0:
                    P_recycled = sum(dest.weight for dest in self.we_driver.recycle_to)
                    self.status_stream.write('%d particles (%g probability) recycled\n' % (n_recycled, P_recycled))
                    log.debug('from: %r' % self.we_driver.recycle_from)
                    log.debug('to:   %r' % self.we_driver.recycle_to)
                    for isrc,src in enumerate(self.we_driver.recycle_from):
                        self.status_stream.write("  {0.count:d} ({0.weight:g} probability) from region {1:d} '{2}'\n"\
                                                 .format(src,isrc,self.system.target_states[isrc].label))
                        
                    for idest, dest in enumerate(self.we_driver.recycle_to):
                        self.status_stream.write(("  {0.count:d} ({0.weight:g} probability,"
                                                  +" {1:%} of recycled particles) to state '{2}'\n")\
                                                 .format(dest, dest.count/n_recycled,
                                                         self.system.initial_states[idest].label))
                else:
                    self.status_stream.write('0 particles recycled\n')
                    
                self.status_stream.write('%d trajectory segments in next iteration\n' % len(next_iter_segments))
                
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
                self.rtracker.end('we_postprocess')
                
                self.rtracker.begin('prep_next_iter')
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
                self.rtracker.end('prep_next_iter')
                
                # Store timing information
                self.rtracker.end('iteration')
                iteration_endtime = self.rtracker.final['iteration'].walltime
                iteration_elapsed = iteration_endtime - iteration_starttime
                                
                iter_walltime = self.rtracker.difference['iteration'].walltime
                iter_cputime = sum(segment.cputime or 0.0 for segment in segments)
                
                iter_summary['walltime'] += iter_walltime
                iter_summary['cputime'] = iter_cputime
                self.data_manager.update_iter_summary(n_iter, iter_summary)
                
                #This may give NaN if starting a truncated simulation
                try:
                    walltime = timedelta(seconds=float(iter_summary['walltime']))
                except ValueError:
                    walltime = 0.0 
                
                try:
                    cputime = timedelta(seconds=float(iter_summary['cputime']))
                except ValueError:
                    cputime = 0.0                     
                    
                self.status_stream.write('Iteration wallclock: {0!s}, cputime: {1!s}\n'\
                                          .format(walltime,
                                                  cputime))
                
    
                # Update HDF5 and flush the status output buffer in preparation for
                # the next iteration
                n_iter += 1
                self.data_manager.current_iteration = n_iter
                self.data_manager.flush_backing()
                self.flush_status()
                
                # Update segments list for next iteration
                segments = next_iter_segments            
                
            # end propagation/WE loop
        finally:
            # Regardless of what happens (end of loop/exit/exception), make sure run data
            # is saved
            self.data_manager.flush_backing()
            self.data_manager.close_backing()

        
        self.status_stream.write('\n%s\n' % time.asctime())
        self.status_stream.write('WEMD run complete.\n')
        self.rtracker.end('run')
        
        # dump resource statistics
        intervals = ResourceUsage(*[timedelta(seconds=field) for field in self.rtracker.difference['run']])
        self.status_stream.write(('\nRun wallclock: {intervals.walltime}, CPU time (WEMD only): {intervals.cputime},'
                                 +' system time: {intervals.systime}\n').format(intervals=intervals))
        
        if wemd.rc.verbose_mode:
            self.status_stream.write('Internal timing information:\n')
            self.rtracker.dump_differences(self.status_stream)
    # end WESimManager.run()
    
    