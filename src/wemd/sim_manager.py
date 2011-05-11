from __future__ import division; __metaclass__ = type

import sys, time, operator, math, numpy
from itertools import izip, imap
from datetime import timedelta
import logging
log = logging.getLogger(__name__)

import wemd
from wemd.util import extloader
from wemd.util.rtracker import ResourceTracker, ResourceUsage
from wemd import Segment
from wemd.util.miscfn import vgetattr

class WESimManager:
    """A state machine and communications broker"""
    def __init__(self, runtime_config = None):
        self.runtime_config = runtime_config or {}
        
        self.data_manager = None
        self.we_driver = None
        self.work_manager = None        
        self.propagator = None
        
        self.system = None
                
        self.status_stream = sys.stdout
        
        self.rtracker = ResourceTracker()                                
            
    def load_data_manager(self):
        drivername = self.runtime_config.get('drivers.data_manager', 'hdf5')
        if drivername.lower() in ('hdf5', 'default'):
            self.data_manager = wemd.data_manager.WEMDDataManager(self)
        else:
            pathinfo = self.runtime_config.get_pathlist('drivers.module_path', default=None)
            self.data_manager = extloader.get_object(drivername, pathinfo)(self)
        log.debug('data manager is %r' % self.data_manager)
            
    def load_we_driver(self):
        drivername = self.runtime_config.get('drivers.we_driver', 'default')
        if drivername.lower() == 'default':
            self.we_driver = wemd.we_driver.WEMDWEDriver(self)
        else:
            pathinfo = self.runtime_config.get_pathlist('drivers.module_path', default=None)
            self.work_manager = extloader.get_object(drivername, pathinfo)(self)
        log.debug('WE algorithm driver is %r' % self.we_driver)
    
    def load_work_manager(self):
        drivername = self.runtime_config.get('args.work_manager_name')
        if not drivername:
            drivername = self.runtime_config.get('drivers.work_manager', 'threads')
        if drivername.lower() == 'serial':
            import wemd.work_managers.serial
            self.work_manager = wemd.work_managers.serial.SerialWorkManager(self)
        elif drivername.lower() == 'processes':
            import wemd.work_managers.processes
            self.work_manager = wemd.work_managers.processes.ProcessWorkManager(self)                    
        elif drivername.lower() in ('threads', 'default'):
            import wemd.work_managers.threads
            self.work_manager = wemd.work_managers.threads.ThreadedWorkManager(self)
        elif drivername.lower() == 'tcpip':
            import wemd.work_managers.tcpip
            self.work_manager = wemd.work_managers.tcpip.TCPWorkManager(self)
        elif drivername.lower() in ('zmq', 'zeromq'):
            import wemd.work_managers.zeromq
            self.work_manager = wemd.work_managers.zeromq.ZMQWorkManager(self)
        else:
            pathinfo = self.runtime_config.get_pathlist('drivers.module_path', default=None)
            self.work_manager = extloader.get_object(drivername, pathinfo)(self)
        log.debug('work manager is %r' % self.work_manager)
        
    def load_propagator(self):
        drivername = self.runtime_config.require('drivers.propagator')
        if drivername.lower() == 'executable':
            import wemd.propagators.executable
            self.propagator = wemd.propagators.executable.ExecutablePropagator(self)
        else:
            pathinfo = self.runtime_config.get_pathlist('drivers.module_path', default=None)
            self.propagator = extloader.get_object(drivername, pathinfo)(self)
        log.debug('propagator is %r' % self.propagator)
        
    def load_system_driver(self):
        sysdrivername = self.runtime_config.require('system.system_driver')
        log.info('loading system driver %r' % sysdrivername)
        pathinfo = self.runtime_config.get_pathlist('system.module_path', default=None)        
        self.system = extloader.get_object(sysdrivername, pathinfo)(self)
        log.debug('system driver is %r' % self.system)
        log.debug('initializing system driver')
        
        # Initialize the system driver
        try:
            self.system.initialize()
        except NotImplementedError:
            msg = 'WEMD systems must now override the initialize() method, not the __init__() method'
            import warnings
            warnings.warn(msg, DeprecationWarning)
            
            # This is really critical, so make sure it gets logged even if DeprecationWarnings are
            # suppressed            
            log.warn(msg)
            
        # Catch a couple loose problems, which will eventually be made illegal but
        # for a transition period will be allowed.
        if self.system.target_states is None:
            msg ='system target_states is None; resetting to []' 
            log.warning(msg)
            warnings.warn(msg, DeprecationWarning)
            self.system.target_states = []

        
    def flush_status(self):
        try:
            self.status_stream.flush()
        except AttributeError:
            pass
        
    def _propagate(self, n_iter, segments):
        # Propagate segments            
        self.rtracker.begin('propagation')
        self.flush_status()
        self.data_manager.close_backing()
        self.work_manager.propagate(segments)
        self.data_manager.open_backing()
        self.data_manager.update_segments(n_iter, segments)
        self.rtracker.end('propagation')
                    
    def run(self):
        """Begin (or continue) running a simulation
        
        If ``mode`` is 'master' (the default), runs the entire WE simulation loop. If
        ``mode`` is 'worker', drops straight into a work manager's "receive-work" loop, if
        supported
        """
        
        if self.work_manager.mode == 'worker':
            self.work_manager.run_worker()
            return
        
        # Set up internal timing
        self.rtracker.begin('run')
        run_starttime = time.time()
        max_walltime = self.runtime_config.get_interval('limits.max_wallclock', default=None, type=float)
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
            max_iter = self.runtime_config.get_int('limits.max_iterations', n_iter+1)
             
            # Guaranteed ordering by seg_id, so segments[seg_id] works for any valid seg_id for this iteration
            segments = self.data_manager.get_segments(n_iter = n_iter)
            while n_iter <= max_iter:
                self.rtracker.begin('iteration')
                if max_walltime and time.time() + 1.1*iteration_elapsed >= run_killtime:
                    self.status_stream.write('Iteration %d would require more than the allotted time. Ending run.\n'
                                             % n_iter)
                    self.work_manager.shutdown(0)
                    sys.exit(0)
                    
                iteration_starttime = time.time()
                self.status_stream.write('\n%s\n' % time.asctime())
                self.status_stream.write('Iteration %d (%d requested)\n' % (n_iter, max_iter))
                
                seg_weights = vgetattr('weight', segments, numpy.float64)
                assert not (seg_weights == 0).any()
                norm = seg_weights.sum()
                self.status_stream.write('norm: %g, error in norm: %g\n' % (norm, norm-1))
                
                segs_to_run = [segment for segment in segments if segment.status != Segment.SEG_STATUS_COMPLETE]
                self.status_stream.write('%d of %d segments remain in iteration %d\n' % (len(segs_to_run), len(segments), n_iter))
                
                if len(segs_to_run) == len(segments):
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
                    n_bins = len(target_counts[target_counts!=0])
                    seg_probs  = vgetattr('weight', segments, numpy.float64)
                    bin_probs  = vgetattr('weight', bins, numpy.float64)
                    min_seg_prob = seg_probs[seg_probs!=0].min()
                    max_seg_prob = seg_probs.max()
                    seg_drange   = math.log(max_seg_prob/min_seg_prob)
                    min_bin_prob = bin_probs[bin_probs!=0].min()
                    max_bin_prob = bin_probs.max()
                    bin_drange = math.log(max_bin_prob/min_bin_prob)
                    n_pop = len(bin_counts[bin_counts!=0])
                    self.status_stream.write('{:d} of {:d} ({:%}) active bins are populated\n'.format(n_pop, n_bins, n_pop/n_bins))
                    self.status_stream.write('per-bin minimum non-zero probability:       {:g}\n'.format(min_bin_prob))
                    self.status_stream.write('per-bin maximum probability:                {:g}\n'.format(max_bin_prob))
                    self.status_stream.write('per-bin probability dynamic range (kT):     {:g}\n'.format(bin_drange))
                    self.status_stream.write('per-segment minimum non-zero probability:   {:g}\n'.format(min_seg_prob))
                    self.status_stream.write('per-segment maximum non-zero probability:   {:g}\n'.format(max_seg_prob))
                    self.status_stream.write('per-segment probability dynamic range (kT): {:g}\n'.format(seg_drange))
                    
                    # Store the above information in HDF5
                    self.data_manager.write_bin_data(n_iter, bin_counts, bin_probs)
                    
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
                    
                    
                    log.debug('preparing work manager for new iteration')
                    self.work_manager.prepare_iteration(n_iter, segments)
                    
                    log.debug('running system-specific per-iteration preprocessing')
                    self.system.preprocess_iteration(n_iter, segments)
                
                # Allow the user to run only one segment to aid in debugging propagators
                log.debug('propagating iteration {:d}'.format(n_iter))
                if self.runtime_config.get('args.only_one_segment',False):
                    log.info('propagating only one segment')
                    self._propagate(n_iter, segs_to_run[0:1])                
                    if len(segs_to_run) > 1:
                        # Exit cleanly after one segment, but only if we haven't finished the last segment in an iteration
                        # In that case, go ahead and run WE
                        break
                else:
                    self._propagate(n_iter, segs_to_run)
                log.debug('propagation complete')
                    
                self.rtracker.begin('we_prep')
                
                # Check to ensure that all segments have been propagated                    
                failed_segments = [segment for segment in segments if segment.status != Segment.SEG_STATUS_COMPLETE]
                if failed_segments:
                    self.status_stream.write('Propagation FAILED for %d segments:\n' % len(failed_segments))
                    for failed_segment in failed_segments:
                        self.status_stream.write('  %d\n' % failed_segment.seg_id)
                    raise RuntimeError('propagation failed for %d segments' % len(failed_segments))
                
                log.debug('running system-specific per-iteration postprocessing')                    
                self.system.postprocess_iteration(n_iter, segments)    
                self.rtracker.end('we_prep')
                
                # Run the weighted ensemble algorithm
                self.rtracker.begin('we_core')
                next_iter_segments = self.we_driver.run_we(segments, self.system.region_set)
                self.rtracker.end('we_core')
    
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
                        log.debug('isrc=%d, src=%r' % (isrc,src))
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
                if self.system.target_states:
                    iter_summary = self.data_manager.get_iter_summary(n_iter)
                    iter_summary['target_flux'] = P_recycled
                    iter_summary['target_hits'] = n_recycled
                    self.data_manager.write_recycling_data(n_iter, self.we_driver.recycle_from)
    
                # Update current segments; their endpoint types were modified by WE                
                self.data_manager.update_segments(n_iter, segments)
                self.data_manager.flush_backing()
                
                log.debug('work manager finalizing iteration')
                self.work_manager.finalize_iteration(n_iter, segments)
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
                        
                # Create new iteration group in HDF5            
                self.data_manager.prepare_iteration(n_iter+1, next_iter_segments, 
                                                    self.system.pcoord_ndim, self.system.pcoord_len, self.system.pcoord_dtype)
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
        
        if self.runtime_config['args.verbose_mode']:
            self.status_stream.write('Internal timing information:\n')
            self.rtracker.dump_differences(self.status_stream)
    # end WESimManager.run()