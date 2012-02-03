"""
HDF5 data manager for WEMD.

Original HDF5 implementation: Joe Kaus
Current implementation: Matt Zwier

WEMD exclusively uses the cross-platform, self-describing file format HDF5
for data storage. This ensures that data is stored efficiently and portably
in a manner that is relatively straightforward for other analysis tools 
(perhaps written in C/C++/Fortran) to access.

"""
from __future__ import division; __metaclass__ = type
import warnings
from operator import attrgetter
from itertools import imap
import numpy
import h5py

import logging
log = logging.getLogger(__name__)

import wemd
from wemd.util.miscfn import vattrgetter
from wemd import Segment

file_format_version = 4

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
                                     ('walltime', numpy.float64)] )

# Using true HDF5 enums here seems to lead to segfaults, so hold off until h5py 
# gets its act together on enums
#seg_status_dtype    = h5py.special_dtype(enum=(numpy.uint8, Segment.statuses))
#seg_initpoint_dtype = h5py.special_dtype(enum=(numpy.uint8, Segment.initpoint_types))
#seg_endpoint_dtype  = h5py.special_dtype(enum=(numpy.uint8, Segment.endpoint_types))
seg_status_dtype = numpy.uint8
seg_initpoint_dtype = numpy.uint8
seg_endpoint_dtype = numpy.uint8


seg_index_dtype = numpy.dtype( [ ('weight', numpy.float64),
                                 ('cputime', numpy.float64),
                                 ('walltime', numpy.float64),
                                 ('parents_offset', numpy.uint32),
                                 ('n_parents', numpy.uint32),                                 
                                 ('status', seg_status_dtype),
                                 ('initpoint_type', seg_initpoint_dtype),
                                 ('endpoint_type', seg_endpoint_dtype), ] )

rec_summary_dtype = numpy.dtype( [ ('count', numpy.uint),
                                   ('weight', numpy.float64) ] )

def warn_deprecated_usage(msg, category=DeprecationWarning, stacklevel=1):
    warnings.warn('deprecated data manager interface use: {}'.format(msg), category, stacklevel+1)

class WEMDDataManager:
    """Data manager for assisiting the reading and writing of WEMD data from/to HDF5 files."""
    
    # defaults for various options
    default_iter_prec = 8
    default_we_h5filename      = 'wemd.h5'
    default_we_h5file_driver   = None
    # Compress any auxiliary dataset whose total size is more than 1MB
    default_aux_compression_threshold = 1048576
        
    def __init__(self):

        #self.h5file = None
        #self.backing_file = backing_file
        #self.system = system
        
        # Backing store filenames; setting these attributes is the only way to specify that
        # a given file should be opened
        self.we_h5filename = self.default_we_h5filename
        self.we_h5file_driver = self.default_we_h5file_driver
        self.h5_access_mode = 'r+'
        self.iter_prec = self.default_iter_prec
        self.aux_compression_threshold = self.default_aux_compression_threshold
        
        # Do not load auxiliary data sets by default, as this can potentially take as much space in RAM
        # as there is auxiliary data stored for a given iteration.
        self.auto_load_auxdata = False
        
        # Backing store h5py.File objects
        self.we_h5file = None
        
        # System description
        self.system = None
        
        # If a configuration file is present, read values from there, otherwise use
        # sensible defaults
        if wemd.rc.config:
            self.we_h5filename  = wemd.rc.config.get_path('data.wemd_data_file', self.default_we_h5filename)
            self.we_h5file_driver = wemd.rc.config.get('data.wemd_data_file_driver', self.default_we_h5file_driver)
            self.iter_prec      = wemd.rc.config.get_int('data.default_iter_prec', self.default_iter_prec)
            self.aux_compression_threshold = wemd.rc.config.get_int('data.default_aux_compression_threshold', 
                                                                    self.default_aux_compression_threshold)
            self.auto_load_auxdata = wemd.rc.config.get_bool('data.load_auxdata', False)
                 
        # A few functions for extracting vectors of attributes from vectors of segments
        self._attrgetters = dict((key, vattrgetter(key)) for key in 
                                 ('seg_id', 'status', 'endpoint_type', 'weight', 'walltime', 'cputime'))
    
    @property
    def h5file(self):
        warn_deprecated_usage('WEMDDataManager.h5file is deprecated in favor of WEMDDataManager.we_h5file')
        return self.we_h5file
        
    def _get_iter_group_name(self, n_iter):
        return 'iter_%0*d' % (self.iter_prec, n_iter)

    def del_iter_group(self, n_iter):
        del self.we_h5file['/iter_%0*d' % (self.iter_prec, n_iter)]

    def get_iter_group(self, n_iter):
        return self.we_h5file['/iter_%0*d' % (self.iter_prec, n_iter)]
        
    @property
    def current_iteration(self):
        return self.we_h5file['/'].attrs['wemd_current_iteration']
    
    @current_iteration.setter
    def current_iteration(self, n_iter):
        self.we_h5file['/'].attrs['wemd_current_iteration'] = n_iter
        
    def open_backing(self):
        '''Open the (already-created) HDF5 file named in self.wemd_h5filename.'''
        if not self.we_h5file:
            # Check to see that the specified file exists and is readable by the HDF5 library
            # throws RuntimeError otherwise; this is not an assert because it should run even
            # when using optimized bytecode (python -O strips "assert" statements).
            h5py.h5f.is_hdf5(self.we_h5filename)
            self.we_h5file = h5py.File(self.we_h5filename, self.h5_access_mode, driver=self.we_h5file_driver)
            try:
                recorded_iter_prec = self.we_h5file['/'].attrs['wemd_iter_prec']
            except KeyError:
                log.info('iteration precision not stored in HDF5; using {:d}'.format(self.iter_prec))
            else:
                log.debug('iteration precision found: {:d}'.format(recorded_iter_prec))
                self.iter_prec = int(recorded_iter_prec)
                    
    def prepare_backing(self):
        '''Create new HDF5 file'''
        
        log.debug('creating HDF5 file {!r} with driver {!r}'.format(self.we_h5filename, self.we_h5file_driver))
        self.we_h5file = h5py.File(self.we_h5filename, 'w', driver=self.we_h5file_driver)
        self.we_h5file['/'].attrs['wemd_file_format_version'] = file_format_version
        self.we_h5file['/'].attrs['wemd_iter_prec'] = self.iter_prec            
        self.current_iteration = 1
        self.we_h5file['/'].create_dataset('summary',
                                        shape=(1,), 
                                        dtype=summary_table_dtype,
                                        maxshape=(None,),
                                        chunks=(100,))
        
    def close_backing(self):
        if self.we_h5file is not None:
            self.we_h5file.close()
            self.we_h5file = None
                    
    def flush_backing(self):
        if self.we_h5file is not None:
            self.we_h5file.flush()
                
    def prepare_iteration(self, n_iter, segments, system=None):
        """Prepare for a new iteration by creating space to store the new iteration's data.
        The number of segments, their IDs, and their lineage must be determined and included
        in the set of segments passed in."""
        
        log.debug('preparing HDF5 group for iteration %d (%d segments)' % (n_iter, len(segments)))
        
        n_particles = len(segments)
        system = system or self.system
        pcoord_ndim = system.pcoord_ndim
        pcoord_len = system.pcoord_len
        pcoord_dtype = system.pcoord_dtype
        n_bins = len(system.region_set.get_all_bins())
        
        # Ensure we have a list for guaranteed ordering
        segments = list(segments)
                
        # Create a table of summary information about each iteration
        summary_table = self.we_h5file['summary']
        if len(summary_table) < n_iter:
            summary_table.resize((n_iter+1,))
        
        iter_group = self.we_h5file.create_group(self._get_iter_group_name(n_iter))
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
        summary_table[n_iter-1] = summary_row
        
        # pcoord is indexed as [particle, time, dimension]
        pcoord_ds = iter_group.create_dataset('pcoord', 
                                              shape=(n_particles, pcoord_len, pcoord_ndim), 
                                              dtype=pcoord_dtype)
        pcoord = pcoord_ds[...]
        
        # Create data sets for tracking trajectories during this iteration
        # Compression is on by default because these are probably rather sparse
        # datasets.
        iter_group.create_dataset('bin_assignments', 
                                  shape=(n_particles, pcoord_len), dtype=numpy.uint32,
                                  chunks=(n_particles, 1), compression='gzip')
        iter_group.create_dataset('bin_populations', 
                                  shape=(pcoord_len, n_bins), dtype=numpy.float64,
                                  chunks=(pcoord_len, 1), compression='gzip')
        iter_group.create_dataset('bin_ntrans',
                                  shape=(n_bins,n_bins),dtype=numpy.uint32,
                                  compression='gzip')
        iter_group.create_dataset('bin_fluxes',
                                  shape=(n_bins,n_bins), dtype=numpy.float64,
                                  compression='gzip')
        iter_group.create_dataset('bin_rates',
                                  shape=(n_bins,n_bins), dtype=numpy.float64,
                                  compression='gzip')
        
        for (seg_id, segment) in enumerate(segments):
            if segment.seg_id is not None:
                assert segment.seg_id == seg_id
            assert segment.p_parent_id is not None
            segment.seg_id = seg_id
            seg_index_table[seg_id]['status'] = segment.status or Segment.SEG_STATUS_UNSET
            seg_index_table[seg_id]['initpoint_type'] = segment.initpoint_type or Segment.SEG_INITPOINT_UNSET
            seg_index_table[seg_id]['endpoint_type']  = segment.endpoint_type or Segment.SEG_ENDPOINT_UNSET
            seg_index_table[seg_id]['weight'] = segment.weight
            seg_index_table[seg_id]['n_parents'] = len(segment.parent_ids)

            # Assign progress coordinate if any exists
            if segment.pcoord is not None:
                if len(segment.pcoord) == 1:
                    # Initial pcoord
                    pcoord[seg_id,0,:] = segment.pcoord[0,:]
                elif segment.pcoord.shape != pcoord.shape[1:]:
                    raise ValueError('segment pcoord shape [%r] does not match expected shape [%r]'
                                     % (segment.pcoord.shape, pcoord.shape[1:]))
                else:
                    pcoord[seg_id,...] = segment.pcoord
                    
        # family tree is stored as two things: a big vector of ints containing parent seg_ids, 
        # and an index (into this vector) and extent pair, stored in the summary table
        
        # voodoo by induction!
        # offset[0] = 0
        # offset[1:] = numpy.add.accumulate(n_parents[:-1])
        seg_index_table[0]['parents_offset'] = 0
        seg_index_table[1:]['parents_offset'] = numpy.add.accumulate(seg_index_table[:-1]['n_parents'])
        
        total_parents = numpy.sum(seg_index_table[:]['n_parents'])
        if total_parents > 0:
            parents_ds = iter_group.create_dataset('parents', (total_parents,), numpy.int32)
            parents = parents_ds[:]
        
            # Don't directly index an HDF5 data set in a loop
            offsets = seg_index_table[:]['parents_offset']
            extents = seg_index_table[:]['n_parents']
            
            for (iseg, segment) in enumerate(segments):
                offset = offsets[iseg]
                extent = extents[iseg]
                assert extent == len(segment.parent_ids)
                assert extent > 0
                assert None not in segment.parent_ids
                assert segment.p_parent_id in segment.parent_ids
                
                # Ensure that the primary parent is first in the list
                parents[offset] = segment.p_parent_id                
                if extent > 1:
                    parent_ids = set(segment.parent_ids)
                    parent_ids.remove(segment.p_parent_id)
                    parent_ids = sorted(parent_ids)
                    if extent == 2:
                        assert len(parent_ids) == 1
                        parents[offset+1] = parent_ids[0]
                    else:
                        parents[offset+1:offset+extent] = parent_ids
                assert set(parents[offset:offset+extent]) == segment.parent_ids
            
            parents_ds[:] = parents                    

        seg_index_table_ds[:] = seg_index_table
        pcoord_ds[...] = pcoord

        del seg_index_table, pcoord
        
    def get_iter_summary(self,n_iter):
        summary_row = numpy.zeros((1,), dtype=summary_table_dtype)
        summary_row[:] = self.we_h5file['summary'][n_iter-1]
        return summary_row
        
    def update_iter_summary(self,n_iter,summary):
        self.we_h5file['summary'][n_iter-1] = summary

    def del_iter_summary(self, min_iter): #delete the iterations starting at min_iter      
        self.we_h5file['summary'].resize((min_iter - 1,))
                     
    def write_recycling_data(self, n_iter, rec_summary):
        iter_group = self.get_iter_group(n_iter)
        rec_data_ds = iter_group.require_dataset('recycling', (len(rec_summary),), dtype=rec_summary_dtype)
        rec_data = numpy.zeros((len(rec_summary),), dtype=rec_summary_dtype)
        for itarget, target in enumerate(rec_summary):
            count, weight = target
            rec_data[itarget]['count'] = count
            rec_data[itarget]['weight'] = weight
        rec_data_ds[...] = rec_data
        
    def write_bin_data(self, n_iter, assignments, populations, n_trans, fluxes, rates):
        iter_group = self.get_iter_group(n_iter)
        try:
            iter_group['bin_assignments'][...] = assignments
            iter_group['bin_populations'][...] = populations
            iter_group['bin_ntrans'][...] = n_trans
            iter_group['bin_fluxes'][...] = fluxes
            iter_group['bin_rates'][...] = rates
        except (KeyError,IndexError,TypeError):
            log.warning('could not write bin data (e.g. assignments/populations/fluxes) for iteration {}'.format(n_iter))
                
    def update_segments(self, n_iter, segments):
        """Update "mutable" fields (status, endpoint type, pcoord, timings, weights) in the HDF5 file
        and update the summary table accordingly.  Note that this DOES NOT update other fields,
        notably the family tree, which is set at iteration preparation and cannot change.
                
        Fields updated:
          * status
          * endpoint type
          * pcoord
          * cputime
          * walltime
          * weight
                    
        Fields not updated:
          * seg_id
          * parents
        """
        
        segments = sorted(list(segments), key=attrgetter('seg_id'))
        seg_ids = [segment.seg_id for segment in segments] 
        
        iter_group = self.get_iter_group(n_iter)
        seg_index_entries = iter_group['seg_index'][seg_ids]
        pcoord_entries = iter_group['pcoord'][seg_ids]
        
        assert len(seg_index_entries) == len(pcoord_entries) == len(seg_ids)
                
        row = numpy.empty((1,), seg_index_dtype)
        for (iseg, segment) in enumerate(segments):
            row[:] = seg_index_entries[iseg]
            row['status'] = segment.status
            row['initpoint_type'] = segment.initpoint_type or Segment.SEG_INITPOINT_UNSET
            row['endpoint_type'] = segment.endpoint_type or Segment.SEG_ENDPOINT_UNSET
            row['cputime'] = segment.cputime
            row['walltime'] = segment.walltime
            row['weight'] = segment.weight
            
            seg_index_entries[iseg] = row
            pcoord_entries[iseg] = segment.pcoord
                    
        iter_group['seg_index'][seg_ids] = seg_index_entries
        iter_group['pcoord'][seg_ids] = pcoord_entries
        
        # Now, to deal with auxiliary data
        # If any segment has any auxiliary data, then the aux dataset must spring into
        # existence. Each is named according to the name in segment.data, and has shape
        # (n_total_segs, ...) where the ... is the shape of the data in segment.data (and may be empty
        # in the case of scalar data) and dtype is taken from the data type of the data entry
        # compression is on by default for datasets that will be more than 1MiB          
        if any(segment.data for segment in segments):
            n_total_segs = len(iter_group['seg_index'])
            aux_group = iter_group.require_group('auxdata')
            for segment in segments:
                for (dsname, data) in segment.data.iteritems():
                    adata = numpy.asarray(data)
                    shape = (n_total_segs,)+adata.shape
                    nbytes = numpy.multiply.reduce(shape)
                    dset = aux_group.require_dataset(dsname,
                                                     shape=shape,
                                                     dtype=adata.dtype,
                                                     compression='gzip' if nbytes > self.aux_compression_threshold else None)
                    dset[segment.seg_id] = adata
                    
            
    def get_segments(self, n_iter, load_auxdata=None):
        '''Return the segments from a given iteration.  This function is optimized for the 
        case of retrieving all segments for a given iteration as quickly as possible, 
        and as such effectively loads all data for the given iteration into memory (which
        is what is currently required for running a WE iteration).
        
        If the optional parameter ``load_auxdata`` is true, then all auxiliary datasets
        available are loaded and mapped onto the ``data`` dictionary of each segment. If 
        ``load_auxdata`` is None, then use the default ``self.auto_load_auxdata``, which can
        be set by the option ``load_auxdata`` in the ``[data]`` section of ``wemd.cfg``. This
        essentially requires as much RAM as there is per-iteration auxiliary data, so this
        behavior is not on by default.'''
        
        if load_auxdata is None: load_auxdata = self.auto_load_auxdata
        
        iter_group = self.get_iter_group(n_iter)
        seg_index_table = iter_group['seg_index'][...]
        pcoords = iter_group['pcoord'][...]
        all_parent_ids = iter_group['parents'][...]
        
        segments = []
        for (seg_id, row) in enumerate(seg_index_table):
            parents_offset = row['parents_offset']
            n_parents = row['n_parents']            
            segment = Segment(seg_id = seg_id,
                              n_iter = n_iter,
                              status = row['status'],
                              n_parents = n_parents,
                              endpoint_type = row['endpoint_type'],
                              walltime = row['walltime'],
                              cputime = row['cputime'],
                              weight = row['weight'],
                              pcoord = pcoords[seg_id])
            parent_ids = all_parent_ids[parents_offset:parents_offset+n_parents]
            segment.p_parent_id = long(parent_ids[0])
            segment.parent_ids = set(imap(long,parent_ids))
            assert len(segment.parent_ids) == n_parents
            segments.append(segment)
        del pcoords
        
        # If any auxiliary data sets are available, load them as well
        if load_auxdata and 'auxdata' in iter_group:
            for (dsname, ds) in iter_group['auxdata'].iteritems():
                for (seg_id, segment) in enumerate(segments):
                    segment.data[dsname] = ds[seg_id]
        
        return segments
    
    def get_segments_by_id(self, n_iter, seg_ids, load_auxdata=None):
        '''Return the given segments from a given iteration.  
        
        If the optional parameter ``load_auxdata`` is true, then all auxiliary datasets
        available are loaded and mapped onto the ``data`` dictionary of each segment. If 
        ``load_auxdata`` is None, then use the default ``self.auto_load_auxdata``, which can
        be set by the option ``load_auxdata`` in the ``[data]`` section of ``wemd.cfg``. This
        essentially requires as much RAM as there is per-iteration auxiliary data, so this
        behavior is not on by default.'''
        
        if load_auxdata is None: load_auxdata = self.auto_load_auxdata

        if len(seg_ids) == 0: return []
        
        iter_group = self.get_iter_group(n_iter)
        seg_index = iter_group['seg_index'][...]
        pcoord_ds = iter_group['pcoord']
        all_parent_ids = iter_group['parents'][...] 
        
        segments = []
        seg_ids = list(seg_ids)
        for seg_id in seg_ids:
            row = seg_index[seg_id]
            parents_offset = row['parents_offset']
            n_parents = row['n_parents']            
            segment = Segment(seg_id = seg_id,
                              n_iter = n_iter,
                              status = row['status'],
                              n_parents = n_parents,
                              endpoint_type = row['endpoint_type'],
                              walltime = row['walltime'],
                              cputime = row['cputime'],
                              weight = row['weight'],)
            parent_ids = all_parent_ids[parents_offset:parents_offset+n_parents]
            segment.p_parent_id = long(parent_ids[0])
            segment.parent_ids = set(imap(long,parent_ids))
            segments.append(segment)
            
        # Use a pointwise selection from pcoord_ds to get only the
        # data we care about
        pcoords_by_seg = pcoord_ds[seg_ids,...]
        for (iseg,segment) in enumerate(segments):
            segment.pcoord = pcoords_by_seg[iseg]
            assert segment.seg_id is not None
            
        if load_auxdata and 'auxdata' in iter_group:
            for (dsname, ds) in iter_group['auxdata'].iteritems():
                for segment in segments:
                    segment.data[dsname] = ds[segment.seg_id]
        
        return segments
        
    def get_children(self, segment):
        '''Return all segments which have the given segment as a parent'''

        if segment.n_iter == self.current_iteration: return []
        
        # Examine the segment index from the following iteration to see who has this segment
        # as a parent.  We don't need to worry about the number of parents each segment
        # has, since each has at least one, and indexing on the offset into the parents array 
        # gives the primary parent ID
        iter_group = self.get_iter_group(segment.n_iter+1)
        all_parent_ids = iter_group['parents'][...]
        seg_index = iter_group['seg_index'][...]
        parent_offsets = seg_index['parents_offset'][...]


        # This is one of the slowest pieces of code I've ever written...
        #seg_index = iter_group['seg_index'][...]
        #seg_ids = [seg_id for (seg_id,row) in enumerate(seg_index) 
        #           if all_parent_ids[row['parents_offset']] == segment.seg_id]
        #return self.get_segments_by_id(segment.n_iter+1, seg_ids)
        p_parents = all_parent_ids[parent_offsets]
        all_seg_ids = numpy.arange(len(parent_offsets), dtype=numpy.uintp)
        seg_ids = all_seg_ids[p_parents == segment.seg_id]
        try:
            len(seg_ids)
        except TypeError:
            seg_ids = [seg_ids]
        
        return self.get_segments_by_id(segment.n_iter+1, seg_ids)

    # The following are dictated by the SimManager interface
    def prepare_run(self):
        self.open_backing()
                
    def finalize_run(self):
        self.close_backing()
        
