"""
HDF5 data manager for WEMD.

Original specifications: Matt Zwier
Original HDF5 implementation: Joe Kaus
Current implementation: Matt Zwier
"""
from __future__ import division; __metaclass__ = type

import operator
import numpy
import h5py

import logging
log = logging.getLogger(__name__)

from wemd.util.miscfn import vattrgetter
from wemd.types import Segment

file_format_version = 1
SUMMARY_TABLE = 'summary'
summary_table_dtype = numpy.dtype( [ ('n_iter', numpy.uint),
                                     ('n_particles', numpy.uint),
                                     ('norm', numpy.float64),
                                     ('target_flux', numpy.float64),
                                     ('target_hits', numpy.uint),
                                     ('min_bin_prob', numpy.float64),
                                     ('max_bin_prob', numpy.float64),
                                     ('bin_dyn_range', numpy.float64),
                                     ('min_seg_prob', numpy.float64),
                                     ('max_seg_prob', numpy.float64),
                                     ('seg_dyn_range', numpy.float64),
                                     ('cputime', numpy.float64),
                                     ('walltime', numpy.float64),
                                     ('status', numpy.byte) ] )
ITER_STATUS_INCOMPLETE = 0
ITER_STATUS_COMPLETE   = 1

seg_index_dtype = numpy.dtype( [ ('weight', numpy.float64),
                                 ('cputime', numpy.float64),
                                 ('walltime', numpy.float64),
                                 ('parents_offset', numpy.uint32),
                                 ('n_parents', numpy.uint32),                                 
                                 ('status', numpy.uint8),
                                 ('endpoint_type', numpy.uint8), ] )
SEG_INDEX_WEIGHT = 0
SEG_INDEX_CPUTIME = 1
SEG_INDEX_WALLTIME = 2
SEG_INDEX_PARENTS_OFFSET = 3
SEG_INDEX_N_PARENTS = 4
SEG_INDEX_STATUS = 5
SEG_INDEX_ENDPOINT_TYPE = 6

rec_summary_dtype = numpy.dtype( [ ('count', numpy.uint),
                                   ('weight', numpy.float64) ] )

class WEMDDataManager:
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
        self.h5file = None
        
        runtime_config = self.sim_manager.runtime_config
        runtime_config.require('data.h5file')
        
        try:
            self.iter_prec = runtime_config.get_int('data_manager.iter_prec')
        except (KeyError,ValueError):
            self.iter_prec = 8 # field width of numeric portion of iteration group names
        
        # A few functions for extracting vectors of attributes from vectors of segments
        self._attrgetters = dict((key, vattrgetter(key)) for key in 
                                 ('seg_id', 'status', 'endpoint_type', 'weight', 'walltime', 'cputime'))
            
    def _get_iter_group_name(self, n_iter):
        return 'iter_%0*d' % (self.iter_prec, n_iter)
    
    def _get_iter_group(self, n_iter):
        return self.h5file['/iter_%0*d' % (self.iter_prec, n_iter)]
    
    @property
    def current_iteration(self):
        return self.h5file['/'].attrs['wemd_current_iteration']
    
    @current_iteration.setter
    def current_iteration(self, n_iter):
        self.h5file['/'].attrs['wemd_current_iteration'] = n_iter
        
    def open_backing(self):
        if not self.h5file:
            self.h5file = h5py.File(self.sim_manager.runtime_config['data.h5file'])
        
    def prepare_backing(self):
        self.open_backing()
        
        self.h5file['/'].attrs['wemd_file_format_version'] = file_format_version
        self.current_iteration = 1
        assert self.h5file['/'].attrs['wemd_current_iteration'] == 1
        
        self.h5file['/'].create_dataset(SUMMARY_TABLE,
                                        shape=(1,), 
                                        dtype=summary_table_dtype,
                                        maxshape=(None,),
                                        chunks=(100,))
        
        
    def close_backing(self):
        self.h5file.close()
        self.h5file = None
        
    def flush_backing(self):
        self.h5file.flush()
                
    def prepare_iteration(self, n_iter, segments, pcoord_ndim, pcoord_len, pcoord_dtype):
        """Prepare for a new iteration by creating space to store the new iteration's data.
        The number of segments, their IDs, and their lineage must be determined and included
        in the set of segments passed in."""
        
        log.debug('preparing HDF5 group for iteration %d (%d segments)' % (n_iter, len(segments)))
        
        n_particles = len(segments)
        
        # Ensure we have a list for guaranteed ordering
        segments = list(segments)
        
        # Create a table of summary information about each iteration
        summary_table = self.h5file[SUMMARY_TABLE]
        if len(summary_table) < n_iter:
            summary_table.resize((n_iter+1,))
        
        iter_group = self.h5file.create_group(self._get_iter_group_name(n_iter))
        iter_group.attrs['n_iter'] = n_iter
        
        # everything indexed by [particle] goes in an index table
        seg_index_table_ds = iter_group.create_dataset('seg_index', shape=(n_particles,),
                                                       dtype=seg_index_dtype)
        # unfortunately, h5py doesn't like in-place modification of individual fields; it expects
        # tuples. So, construct everything in a numpy array and then dump the whole thing into hdf5
        # In fact, this appears to be an h5py best practice (collect as much in ram as possible and then dump)
        seg_index_table = numpy.zeros((n_particles,), dtype=seg_index_dtype)
                
        summary_row = numpy.zeros((1,), dtype=summary_table_dtype)
        summary_row['n_iter'] = n_iter
        summary_row['n_particles'] = n_particles
        summary_row['norm'] = numpy.add.reduce(map(self._attrgetters['weight'], segments))
        summary_row['status'] = ITER_STATUS_INCOMPLETE
        summary_table[n_iter-1] = summary_row
        
        # pcoord is indexed as [particle, time, dimension]
        pcoord_ds = iter_group.create_dataset('pcoord', 
                                              shape=(n_particles, pcoord_len, pcoord_ndim), 
                                              dtype=pcoord_dtype)
        pcoord = pcoord_ds[...]
        log.debug('pcoord shape is {!r}'.format(pcoord.shape))
        
        for (seg_id, segment) in enumerate(segments):
            if log.isEnabledFor(logging.DEBUG):
                log.debug('processing segment %r' % segment)
                log.debug('assigning seg_id=%r' % seg_id)
            assert segment.p_parent_id is not None
            segment.seg_id = seg_id
            seg_index_table[seg_id]['status'] = segment.status
            seg_index_table[seg_id]['weight'] = segment.weight
            seg_index_table[seg_id]['n_parents'] = len(segment.parent_ids) if segment.parent_ids else 1

            # Assign progress coordinate if any exists
            if segment.pcoord is not None:
                if len(segment.pcoord) == 1:
                    # Initial pcoord
                    pcoord[seg_id,0,:] = segment.pcoord[0,:]
                else:
                    raise ValueError('segment pcoord shape [%r] does not match expected shape [%r]'
                                     % (segment.pcoord.shape, pcoord.shape[1:]))
            

                
        # family tree is stored as two things: a big vector of ints containing parent seg_ids, 
        # and an index (into this vector) and extent pair
        
        # voodoo by induction!
        # offset[0] = 0
        # offset[1:] = numpy.add.accumulate(n_parents[:-1])
        seg_index_table[0]['parents_offset'] = 0
        seg_index_table[1:]['parents_offset'] = numpy.add.accumulate(seg_index_table[:-1]['n_parents'])
        
        total_parents = numpy.sum(seg_index_table[:]['n_parents'])
        if total_parents > 0:
            log.debug('creating dataset for %d parents' % total_parents)
            parents_ds = iter_group.create_dataset('parents', (total_parents,), numpy.int32)
            parents = parents_ds[:]
        
            # Don't directly index an HDF5 data set in a loop
            offsets = seg_index_table[:]['parents_offset']
            extents = seg_index_table[:]['n_parents']
            
            for (iseg, segment) in enumerate(segments):
                offset = offsets[iseg]
                extent = extents[iseg]
                assert extent == len(segment.parent_ids)
                #log.debug('segment %d has %d parent(s): %r'  % (iseg, extent, tuple(sorted(segment.parent_ids))))
                if extent == 0: continue
                
                # Ensure that the primary parent is first in the list
                # If any parent is None, indicating that a recycled particle was merged in and selected as
                # the initial position of a segment, store -1 to indicate this
                assert segment.p_parent_id in segment.parent_ids
                parents[offset] = segment.p_parent_id if segment.p_parent_id is not None else -1
                segment.parent_ids.remove(segment.p_parent_id)
                if None in segment.parent_ids:
                    segment.parent_ids.remove(None)
                    segment.parent_ids.add(-1)
                if extent == 2:
                    parents[offset+1] = segment.parent_ids.pop()
                else: # extent > 2:
                    parents[offset+1:offset+extent] = list(segment.parent_ids)
            
            parents_ds[:] = parents                    

        # Since we accumulated many of these changes in RAM (and not directly in HDF5), propagate
        # the changes out to HDF5
        seg_index_table_ds[:] = seg_index_table

        pcoord_ds[...] = pcoord
        #self.flush_backing()
    
    def update_segments(self, n_iter, segments):
        """Update "mutable" fields (status, endpoint type, pcoord, timings) in the HDF5 file
        and update the summary table accordingly.  Note that this DOES NOT update other fields,
        notably weights and the family tree.
                
        Fields updated:
          * status
          * endpoint type
          * pcoord
          * cputime
          * walltime
          
        Fields not updated:
          * seg_id
          * parents
          * weight
        """
        
        segments = list(segments)
        iter_group = self._get_iter_group(n_iter)
        seg_index_table = iter_group['seg_index'][...]
        pcoords = iter_group['pcoord'][...]
        
        row = numpy.empty((1,), seg_index_dtype)
        for segment in segments:
            seg_id = segment.seg_id
            row[:] = seg_index_table[seg_id]
            #log.debug('row: %r' % (row,))
            row['status'] = segment.status
            row['endpoint_type'] = segment.endpoint_type or Segment.SEG_ENDPOINT_TYPE_NOTSET
            row['cputime'] = segment.cputime
            row['walltime'] = segment.walltime
            #log.debug('row: %r' % (row,))
            seg_index_table[seg_id] = row
            
            pcoords[seg_id] = segment.pcoord
        
        iter_group['seg_index'][...] = seg_index_table
        iter_group['pcoord'][...] = pcoords
        
    def load_parents(self, segments):
        """Load the parents of the given segments, returning a mapping of seg_id => (p_parent, set(parents)) """
        raise NotImplementedError
        
    def get_segments(self, n_iter, status = None, endpoint_type = None):
        iter_group = self._get_iter_group(n_iter)
        seg_index_table = iter_group['seg_index'][...]
        
        segments = []
        for (seg_id, row) in enumerate(seg_index_table):
            if status is not None and row['status'] != status: continue
            if endpoint_type is not None and row['endpoint_type'] != endpoint_type: continue
            
            segment = Segment(seg_id = seg_id,
                              n_iter = n_iter,
                              status = row['status'],
                              n_parents = row['n_parents'],
                              endpoint_type = row['endpoint_type'],
                              walltime = row['walltime'],
                              cputime = row['cputime'],
                              weight = row['weight'])
            
            segment.pcoord = iter_group['pcoord'][seg_id]
            segments.append(segment)
        return segments
    
    def get_iter_summary(self,n_iter):
        summary_row = numpy.zeros((1,), dtype=summary_table_dtype)
        summary_row[:] = self.h5file[SUMMARY_TABLE][n_iter-1]
        return summary_row
        
    def update_iter_summary(self,n_iter,summary):
        self.h5file[SUMMARY_TABLE][n_iter-1] = summary
        
    def write_bin_data(self, n_iter, bin_counts, bin_probabilities):
        assert len(bin_counts) == len(bin_probabilities)
        iter_group = self._get_iter_group(n_iter)
        bin_count_ds = iter_group.require_dataset('bin_counts', (len(bin_counts),), numpy.uint)
        bin_prob_ds  = iter_group.require_dataset('bin_probs', (len(bin_probabilities),), numpy.float64)
        
        bin_count_ds[:] = bin_counts
        bin_prob_ds[:]  = bin_probabilities
        
    def write_recycling_data(self, n_iter, rec_summary):
        iter_group = self._get_iter_group(n_iter)
        rec_data_ds = iter_group.require_dataset('recycling', (len(rec_summary),), dtype=rec_summary_dtype)
        rec_data = numpy.zeros((len(rec_summary),), dtype=rec_summary_dtype)
        for itarget, target in enumerate(rec_summary):
            count, weight = target
            rec_data[itarget]['count'] = count
            rec_data[itarget]['weight'] = weight
        rec_data_ds[...] = rec_data
        

        
            
        
        