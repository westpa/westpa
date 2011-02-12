from __future__ import division; __metaclass__ = type

import multiprocessing, logging, math, itertools, numpy

log = logging.getLogger(__name__)

from wemd.work_managers import WEMDWorkManager

class ProcessWorkManager(WEMDWorkManager):
    def __init__(self, sim_manager):
        log.debug('initializing threaded work manager')
        super(ProcessWorkManager,self).__init__(sim_manager)
        self.cpu_count = multiprocessing.cpu_count()
        log.debug('cpu count: %d' % self.cpu_count)
        self.n_procs = sim_manager.runtime_config.get_int('work_manager.n_threads', self.cpu_count)
        self.n_iter = None
    
        log.info('using %d processes for parallel propagation' % self.n_procs)
        
    def propagate(self, segments):
        # Determine the shape of a square array with n_procs rows large enough
        # to hold all our segments.  After assigning segments to this array,
        # slicing along rows can be used to assign sets of segments to workers
        n_rounds = int(math.ceil(len(segments) / self.n_procs))
        segarray = numpy.empty((self.n_procs, n_rounds), numpy.object_)
        segarray.flat[0:len(segments)] = segments
        
        # Segments are assigned to processes here, and they are carried across
        # the call to fork().  Queues are used to return segments to the
        # master.        
        queues = [multiprocessing.Queue() for i in xrange(0, self.n_procs)]
        processes = [WorkerProcess(self.sim_manager, segarray[iproc,:], queues[iproc])
                     for iproc in xrange(0, self.n_procs)]
            
        for process in processes:
            # Spawn processes, begin propagation in each
            process.start()
            
        # Might as well do a little work here
        segs_by_id = {segment.seg_id: segment for segment in segments}
        
        # Wait on all processes
        for process in processes: 
            process.join()      
        
        # Get all the segments from the queue
        propagated_segs = list(itertools.chain(*[queue.get() for queue in queues]))
        
        # Check to make sure we got back what we sent out
        out_ids = set(segment.seg_id for segment in segments)
        in_ids  = set(segment.seg_id for segment in propagated_segs if segment is not None)
        assert out_ids == in_ids
        
        for propagated_seg in propagated_segs:
            orig_segment = segs_by_id[propagated_seg.seg_id]
            orig_segment.status = propagated_seg.status
            orig_segment.walltime = propagated_seg.walltime
            orig_segment.cputime  = propagated_seg.cputime
            orig_segment.pcoord[...] = propagated_seg.pcoord[...]
        
class WorkerProcess(multiprocessing.Process):
    def __init__(self, sim_manager, segments, queue):
        multiprocessing.Process.__init__(self)
        self.sim_manager = sim_manager
        self.segments = segments
        self.queue = queue
        
    def run(self):
        propagator = self.sim_manager.propagator
        
        log.debug('propagating %d segments' % len(self.segments))
        for i in xrange(0,len(self.segments)):
            segment = self.segments[i]
   
            if segment is None: continue

            propagator.propagate([segment])
            
            self.segments[i] = segment
            
        self.queue.put(self.segments)
        log.debug('propagation complete')
