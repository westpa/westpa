"""
HDF5 data manager for WEMD.

Original HDF5 implementation: Joe Kaus
Current implementation: Matt Zwier

WEMD exclusively uses the cross-platform, self-describing file format HDF5
for data storage. This ensures that data is stored efficiently and portably
in a manner that is relatively straightforward for other analysis tools 
(perhaps written in C/C++/Fortran) to access.

The data is laid out in HDF5 as follows:
    - summary -- overall summary data for the simulation
    - /iterations/ -- data for individual iterations, one group per iteration under /iterations
        - iter_00000001/ -- data for iteration 1
            - seg_index -- overall information about segments in the iteration, including weight
            - pcoord -- progress coordinate data organized as [seg_id][time][dimension]
            - parents -- data used to reconstruct the split/merge history of trajectories
            - recycling -- flux and event count for recycled particles, on a per-target-state basis
            - aux_data/ -- auxiliary datasets (data stored on the 'data' field of Segment objects)
    - /basis_states -- data for basis states
        - index -- index to basis states (weights, names, when they are valid in the simulation history, etc)
        - pcoord_00000000 -- progress coordinate data for first set of basis states
    - /initial_states -- data for initial/recycling states, generated from basis states
        - index -- index to basis states (weights, names, when they were generated/used, etc); note that states for
                   which iteration_generated == 0 are used for initial points, and those states with "0" for iteration_used
                   have not been used.
        - pcoord -- progress coordinate data for initial/recycling states
    - /target_states -- data for target states
        - index -- index to target states (names and when they are valid in the simulation history)
        - pcoords -- progress coordinates which map to the bins of target states

The file root object has an integer attribute 'wemd_file_format_version' which can be used to
determine how to access data even as the file format (i.e. organization of data within HDF5 file)
evolves. 

Version history:
    Version 5 
        - moved iter_* groups into a top-level iterations/ group,
        - added in-HDF5 storage for basis states, target states, and generated states
"""        
from __future__ import division; __metaclass__ = type
import warnings
from operator import attrgetter
from itertools import imap, izip
import numpy
import h5py
import threading

import logging
log = logging.getLogger(__name__)

import wemd
from wemd.util.miscfn import vattrgetter
from wemd.segment import Segment
from wemd.states import BasisState, TargetState

file_format_version = 5

class dummy_lock:
    def __init__(self):
        self.depth = 0
        self.tid = threading.current_thread().ident
    
    def __enter__(self):
        assert threading.current_thread().ident == self.tid
        self.depth += 1
        log.debug('acquired data manager lock (depth now {:d})'.format(self.depth))
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        assert threading.current_thread().ident == self.tid
        self.depth -= 1
        log.debug('returned data manager lock (depth now {:d})'.format(self.depth))
        
class flushing_lock:
    def __init__(self, lock, fileobj):
        self.lock = lock
        self.fileobj = fileobj
        
    def __enter__(self):
        self.lock.acquire()
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fileobj.flush()
        self.lock.release()

# Data types for use in the HDF5 file

vstr_dtype = h5py.new_vlen(str)
h5ref_dtype = h5py.special_dtype(ref=h5py.Reference)

# Using true HDF5 enums here seems to lead to segfaults, so hold off until h5py 
# gets its act together on enums
#seg_status_dtype    = h5py.special_dtype(enum=(numpy.uint8, Segment.statuses))
#seg_initpoint_dtype = h5py.special_dtype(enum=(numpy.uint8, Segment.initpoint_types))
#seg_endpoint_dtype  = h5py.special_dtype(enum=(numpy.uint8, Segment.endpoint_types))
seg_status_dtype = numpy.uint8
seg_initpoint_dtype = numpy.uint8
seg_endpoint_dtype = numpy.uint8
    
summary_table_dtype = numpy.dtype( [ ('n_particles', numpy.uint),      # Number of live trajectories in this iteration
                                     ('norm', numpy.float64),          # Norm of probability, to watch for errors or drift
                                     ('min_bin_prob', numpy.float64),  # Per-bin minimum probability
                                     ('max_bin_prob', numpy.float64),  # Per-bin maximum probability
                                     ('min_seg_prob', numpy.float64),  # Per-segment minimum probability
                                     ('max_seg_prob', numpy.float64),  # Per-segment maximum probability
                                     ('cputime', numpy.float64),       # Total CPU time for this iteration
                                     ('walltime', numpy.float64)] )    # Total wallclock time for this iteration


seg_index_dtype = numpy.dtype( [ ('weight', numpy.float64),               # Statistical weight of this segment
                                 ('cputime', numpy.float64),              # CPU time used in propagating this segment
                                 ('walltime', numpy.float64),             # Wallclock time used in propagating this segment
                                 ('parents_offset', numpy.uint32),  # Offset into the parents array for the parents of this segment
                                 ('n_parents', numpy.uint32),             # Number of parents this segment has
                                 ('status', seg_status_dtype),            # Status of propagation of this segment
                                 ('initpoint_type', seg_initpoint_dtype), # Initial point type (basis/initial state or continuation)
                                 ('endpoint_type', seg_endpoint_dtype),   # Endpoint type (will continue, merged, or recycled)
                                 ] )

rec_summary_dtype = numpy.dtype( [ ('count', numpy.uint),        # Number of transitions into this state
                                   ('flux', numpy.float64) ] )   # Probability flux into this state


bstate_index_dtype = numpy.dtype( [('iter_created', numpy.uint),# When this basis state set was created (i.e. when it is valid)
                                   ('group_ref', h5ref_dtype)]) # A reference to the group containing the state info

bstate_dtype = numpy.dtype( [ ('label', vstr_dtype),            # An optional descriptive label
                              ('probability', numpy.float64),   # Probability that this state will be selected 
                              ('auxref', vstr_dtype),           # An optional auxiliar data reference 
                              ])


# Even when initial state generation is off and basis states are passed through directly, an initial state entry
# is created, as that allows precise tracing of where the initial state information came from (since basis states can
# change from iteration to iteration, only the combination of iter_created and basis_state_id is unique, and storing
# that here avoids having to put more fields in the segment index to track that case).

istate_dtype = numpy.dtype( [('iter_created', numpy.uint),      # Iteration during which this state was generated (0 for at w_init)
                             ('iter_used', numpy.uint),         # When this state was used to start a new trajectory 
                                                                # (1 for iteration 1, etc)
                             ('basis_state_id', numpy.uint)])   # Which basis state this state was generated from

tstate_index_dtype = numpy.dtype([('iter_created', numpy.uint), # Iteration when this state list is valid
                                  ('group_ref', h5ref_dtype)])  # Reference to a group containing further data; this will be the
                                                                # null reference if there is no target state for that timeframe.
tstate_dtype = numpy.dtype( [('label', vstr_dtype),])           # An optional descriptive label for this state
                             
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
    
    def flushing_lock(self):
        return flushing_lock(self.lock, self.we_h5file)
        
    def __init__(self):        
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
        
        # If a configuration file is present, read values from there, otherwise use
        # sensible defaults
        if wemd.rc.config:
            self.we_h5filename  = wemd.rc.config.get_path('data.wemd_data_file', self.default_we_h5filename)
            self.we_h5file_driver = wemd.rc.config.get('data.wemd_data_file_driver', self.default_we_h5file_driver)
            self.iter_prec      = wemd.rc.config.get_int('data.default_iter_prec', self.default_iter_prec)
            self.aux_compression_threshold = wemd.rc.config.get_int('data.default_aux_compression_threshold', 
                                                                    self.default_aux_compression_threshold)
            self.auto_load_auxdata = wemd.rc.config.get_bool('data.load_auxdata', False)
        self.lock = threading.RLock()
                 
        # A few functions for extracting vectors of attributes from vectors of segments
        self._attrgetters = dict((key, vattrgetter(key)) for key in 
                                 ('seg_id', 'status', 'endpoint_type', 'weight', 'walltime', 'cputime'))
    
    @property
    def h5file(self):
        warn_deprecated_usage('WEMDDataManager.h5file is deprecated in favor of WEMDDataManager.we_h5file')
        return self.we_h5file

    def iter_group_name(self, n_iter, absolute=True):
        if absolute:
            return '/iterations/iter_{:0{prec}d}'.format(n_iter, prec=self.iter_prec)
        else:
            return 'iter_{:0{prec}d}'.format(n_iter, prec=self.iter_prec)

    def create_iter_group(self, n_iter):
        with self.lock:
            iter_group = self.we_h5file.create_group('/iterations/iter_{:0{prec}d}'.format(n_iter, prec=self.iter_prec))
            iter_group.attrs['n_iter'] = n_iter
        return iter_group
            
    def del_iter_group(self, n_iter):
        with self.lock:
            del self.we_h5file['/iterations/iter_{:0{prec}d}'.format(n_iter, prec=self.iter_prec)]

    def get_iter_group(self, n_iter):
        with self.lock:
            return self.we_h5file['/iterations/iter_{:0{prec}d}'.format(n_iter, prec=self.iter_prec)]
            
    def get_seg_index(self, n_iter):
        with self.lock:
            seg_index = self.get_iter_group(n_iter)['seg_index']
            return seg_index
        
    @property
    def current_iteration(self):
        with self.lock:
            return self.we_h5file['/'].attrs['wemd_current_iteration']
    
    @current_iteration.setter
    def current_iteration(self, n_iter):
        with self.lock:
            self.we_h5file['/'].attrs['wemd_current_iteration'] = n_iter
        
    def open_backing(self):
        '''Open the (already-created) HDF5 file named in self.wemd_h5filename.'''
        if not self.we_h5file:
            # Check to see that the specified file exists and is readable by the HDF5 library
            # throws RuntimeError otherwise; this is not an assert because it should run even
            # when using optimized bytecode (python -O strips "assert" statements).
            self.we_h5file = h5py.File(self.we_h5filename, self.h5_access_mode, driver=self.we_h5file_driver)
            try:
                recorded_iter_prec = self.we_h5file['/'].attrs['wemd_iter_prec']
            except KeyError:
                log.info('iteration precision not stored in HDF5; using {:d}'.format(self.iter_prec))
            else:
                log.debug('iteration precision found: {:d}'.format(recorded_iter_prec))
                self.iter_prec = int(recorded_iter_prec)
                    
    def prepare_backing(self): #istates):
        '''Create new HDF5 file'''
        self.we_h5file = h5py.File(self.we_h5filename, 'a', driver=self.we_h5file_driver)
        
        with self.flushing_lock():
            self.we_h5file['/'].attrs['wemd_file_format_version'] = file_format_version
            self.current_iteration = 0
            self.we_h5file['/'].create_dataset('summary',
                                               shape=(1,), 
                                               dtype=summary_table_dtype,
                                               maxshape=(None,))
            self.we_h5file.create_group('/iterations')
        
    def close_backing(self):
        if self.we_h5file is not None:
            with self.lock:
                self.we_h5file.close()
            self.we_h5file = None
        
    def flush_backing(self):
        if self.we_h5file is not None:
            with self.lock:
                self.we_h5file.flush()

    def save_basis_states(self, bstates):
        '''Save the given basis states in the HDF5 file; they will be used for the next iteration to
        be propagated.  A complete set is required, even if nominally appending to an existing set;
        no processing is done of previously-existing basis states except to update their 
        valid-until field, which simplifies the mapping of IDs to the table.'''
        
        system = wemd.rc.get_system_driver()
        
        # Assemble all the important data before we start to modify the HDF5 file
        bstates = list(bstates)
        assert bstates, 'non-empty state list required'
        state_table = numpy.empty((len(bstates),), dtype=bstate_dtype)
        state_pcoords = numpy.empty((len(bstates),system.pcoord_ndim), dtype=system.pcoord_dtype)
        for i, state in enumerate(bstates):
            state.state_id = i
            state_table[i]['label'] = state.label
            state_table[i]['probability'] = state.probability
            state_table[i]['auxref'] = state.auxref or ''
            state_pcoords[i] = state.pcoord
        
        # Commit changes to HDF5
        with self.lock:
            try:
                master_bstates_group = self.we_h5file['/basis_states']
            except KeyError:
                master_bstates_group = self.we_h5file.create_group('/basis_states')
                master_index = master_bstates_group.create_dataset('index', shape=(1,), maxshape=(None,), dtype=bstate_index_dtype)
                n_sets = 1
            else:
                master_index = master_bstates_group['index'][...]
                n_sets = len(master_index) + 1
                master_index.resize((n_sets,))
            
            n_iter = self.current_iteration
            state_group = master_bstates_group.create_group(self.iter_group_name(n_iter, absolute=False))                
            master_index_row = numpy.empty((), dtype=bstate_index_dtype)
            master_index_row['iter_created'] = n_iter
            master_index_row['group_ref'] = state_group.ref            
            state_group['index'] = state_table
            state_group['pcoord'] = state_pcoords
            master_bstates_group['index'][n_sets-1] = master_index_row

    def get_basis_states(self, n_iter):
        '''Return a list of BasisState objects representing the basis states that are in use for iteration n_iter.'''
        
        with self.lock:
            master_index = self.we_h5file['/basis_states/index'][...]
            
            i = 0
            while i < len(master_index) and master_index[i]['iter_created'] < n_iter:
                i += 1
            
            bstate_group = self.we_h5file[master_index[i]['group_ref']]
            bstate_index = bstate_group['index']
            bstate_pcoords = bstate_group['pcoord']
            
            bstates = [BasisState(state_id=i, label=row['label'], probability=row['probability'],
                                  auxref = row['auxref'] or None, pcoord=pcoord)
                       for (i, (row, pcoord))  in enumerate(izip(bstate_index, bstate_pcoords))]
            return bstates
            
    def save_target_states(self, tstates):
        '''Save the given basis states in the HDF5 file; they will be used for the next iteration to
        be propagated.  A complete set is required, even if nominally appending to an existing set;
        no processing is done of previously-existing basis states except to update their 
        valid-until field, which simplifies the mapping of IDs to the table.'''
        
        system = wemd.rc.get_system_driver()

        # Assemble all the important data before we start to modify the HDF5 file
        tstates = list(tstates)
        if tstates:            
            state_table = numpy.empty((len(tstates),), dtype=tstate_dtype)
            state_pcoords = numpy.empty((len(tstates),system.pcoord_ndim), dtype=system.pcoord_dtype)
            for i, state in enumerate(tstates):
                state.state_id = i
                state_table[i]['label'] = state.label
                state_pcoords[i] = state.pcoord
        else:
            state_table = None
            state_pcoords = None
        
        # Commit changes to HDF5
        with self.lock:
            try:
                master_tstates_group = self.we_h5file['/target_states']
            except KeyError:
                master_tstates_group = self.we_h5file.create_group('/target_states')
                master_index = master_tstates_group.create_dataset('index', shape=(1,), maxshape=(None,), dtype=tstate_index_dtype)
                n_sets = 1
            else:
                master_index = master_tstates_group['index'][...]
                n_sets = len(master_index) + 1
                master_index.resize((n_sets,))
            
            n_iter = self.current_iteration
            master_index_row = numpy.empty((), dtype=bstate_index_dtype)
            master_index_row['iter_created'] = n_iter
            
            if tstates:
                state_group = master_tstates_group.create_group(self.iter_group_name(n_iter, absolute=False))                
                master_index_row['group_ref'] = state_group.ref            
                state_group['index'] = state_table
                state_group['pcoord'] = state_pcoords
            else:
                master_index_row['group_ref'] = None # this is the sentinel that no target states exist for this time frame
                
            master_tstates_group['index'][n_sets-1] = master_index_row

    def get_target_states(self, n_iter):
        '''Return a list of Target objects representing the target (sink) states that are in use for iteration n_iter.'''
        
        with self.lock:
            master_index = self.we_h5file['/target_states/index'][...]
            
            i = 0
            while i < len(master_index) and master_index[i]['iter_created'] < n_iter:
                i += 1
            
            group_ref = master_index[i]['group_ref']
            if not group_ref:
                # No states for this iteration
                return []
            else:
                tstate_group = self.we_h5file[master_index[i]['group_ref']]
                tstate_index = tstate_group['index']
                tstate_pcoords = tstate_group['pcoord']
                
                tstates = [TargetState(state_id=i, label=row['label'], pcoord=pcoord)
                           for (i, (row, pcoord))  in enumerate(izip(tstate_index, tstate_pcoords))]
                return tstates
            
    def alloc_initial_states(self,n_states):
        '''Reserve space for ``n_states`` initial states, and return their IDs'''
        
        system = wemd.rc.get_system_driver()
        with self.lock:
            try:
                istate_group = self.we_h5file['/initial_states']
            except KeyError:
                istate_group = self.we_h5file.create_group('/initial_states')
                index_ds = istate_group.create_dataset('index', shape=(n_states,), maxshape=(None,),
                                                       dtype=istate_dtype)
                pcoord_ds = istate_group.create_dataset('pcoord', 
                                                        shape=(n_states,system.pcoord_ndim),
                                                        maxshape=(None,system.pcoord_ndim),
                                                        dtype=system.pcoord_dtype)
                n_current_states = 0
            else:
                index_ds = istate_group['index']
                pcoord_ds = istate_group['pcoord']
                n_current_states = len(index_ds)
                index_ds.resize((n_current_states+n_states,))
                pcoord_ds.resize((n_current_states+n_states,system.pcoord_ndim))
        return range(n_current_states,n_current_states+n_states)
            
    def save_initial_states(self, initial_states):
        '''Save the given initial states in the HDF5 file'''
        
        system = wemd.rc.get_system_driver()
        
        initial_states = sorted(initial_states,key=attrgetter('state_id'))
        state_ids = [state.state_id for state in initial_states]
        index_entries = numpy.empty((len(initial_states),), dtype=istate_dtype)
        pcoord_vals = numpy.empty((len(initial_states), system.pcoord_ndim), dtype=system.pcoord_dtype)
        for i, initial_state in enumerate(initial_states):
            index_entries[i]['iter_created'] = initial_state.iter_created
            index_entries[i]['iter_used'] = initial_state.iter_used or 0
            index_entries[i]['basis_state_id'] = initial_state.basis_state_id
            pcoord_vals[i] = initial_state.pcoord
        
        with self.lock:
            init_state_group = self.we_h5file['/initial_states']
            init_state_group['index'][state_ids] = index_entries
            init_state_group['pcoord'][state_ids] = pcoord_vals        
                
    def prepare_iteration(self, n_iter, segments):
        """Prepare for a new iteration by creating space to store the new iteration's data.
        The number of segments, their IDs, and their lineage must be determined and included
        in the set of segments passed in."""
        
        log.debug('preparing HDF5 group for iteration %d (%d segments)' % (n_iter, len(segments)))
        
        # Ensure we have a list for guaranteed ordering
        segments = list(segments)
        n_particles = len(segments)
        system = wemd.rc.get_system_driver()
        pcoord_ndim = system.pcoord_ndim
        pcoord_len = system.pcoord_len
        pcoord_dtype = system.pcoord_dtype
        n_bins = len(system.new_region_set().get_all_bins())
        
        with self.lock:
            # Create a table of summary information about each iteration
            summary_table = self.we_h5file['summary']
            if len(summary_table) < n_iter:
                summary_table.resize((n_iter+1,))
            
            iter_group = self.create_iter_group(n_iter)
            iter_group.attrs['n_iter'] = n_iter
            
            # everything indexed by [particle] goes in an index table
            seg_index_table_ds = iter_group.create_dataset('seg_index', shape=(n_particles,),
                                                           dtype=seg_index_dtype)
            # unfortunately, h5py doesn't like in-place modification of individual fields; it expects
            # tuples. So, construct everything in a numpy array and then dump the whole thing into hdf5
            # In fact, this appears to be an h5py best practice (collect as much in ram as possible and then dump)
            seg_index_table = numpy.zeros((n_particles,), dtype=seg_index_dtype)
                    
            summary_row = numpy.zeros((1,), dtype=summary_table_dtype)
            summary_row['n_particles'] = n_particles
            summary_row['norm'] = numpy.add.reduce(map(self._attrgetters['weight'], segments))
            summary_table[n_iter-1] = summary_row
            
            # pcoord is indexed as [particle, time, dimension]
            pcoord_ds = iter_group.create_dataset('pcoord', 
                                                  shape=(n_particles, pcoord_len, pcoord_ndim), 
                                                  dtype=pcoord_dtype)
            pcoord = pcoord_ds[...]
            
            _assignments_ds = iter_group.create_dataset('bin_assignments', shape=(n_particles, pcoord_len), dtype=numpy.uint32)
            _populations_ds = iter_group.create_dataset('bin_populations', shape=(pcoord_len, n_bins), dtype=numpy.float64)
            _n_trans_ds = iter_group.create_dataset('bin_ntrans', shape=(n_bins,n_bins), dtype=numpy.uint32)
            _fluxes_ds = iter_group.create_dataset('bin_fluxes', shape=(n_bins,n_bins), dtype=numpy.float64)
            _rates_ds = iter_group.create_dataset('bin_rates', shape=(n_bins,n_bins), dtype=numpy.float64)
            
            for (seg_id, segment) in enumerate(segments):
                if segment.seg_id is not None:
                    assert segment.seg_id == seg_id
                else:
                    segment.seg_id = seg_id
                # Parent must be set, though what it means depends on initpoint_type
                assert segment.p_parent_id is not None
                assert segment.initpoint_type is not None
                segment.seg_id = seg_id
                seg_index_table[seg_id]['status'] = segment.status
                seg_index_table[seg_id]['weight'] = segment.weight
                seg_index_table[seg_id]['n_parents'] = len(segment.parent_ids)
                seg_index_table[seg_id]['initpoint_type'] = segment.initpoint_type
    
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
            # and an index (into this vector) and extent pair
            
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
                        parent_ids = list(sorted(parent_ids))
                        if extent == 2:
                            assert len(parent_ids) == 1
                            parents[offset+1] = parent_ids[0]
                        else:
                            parents[offset+1:offset+extent] = parent_ids
                    assert set(parents[offset:offset+extent]) == segment.parent_ids
                
                parents_ds[:] = parents                    
    
            # Since we accumulated many of these changes in RAM (and not directly in HDF5), propagate
            # the changes out to HDF5
            seg_index_table_ds[:] = seg_index_table
            pcoord_ds[...] = pcoord
            
            # A few explicit deletes
            del seg_index_table, pcoord
        
    
    def get_iter_summary(self,n_iter):
        with self.lock:
            summary_row = numpy.zeros((1,), dtype=summary_table_dtype)
            summary_row[:] = self.we_h5file['summary'][n_iter-1]
            return summary_row
        
    def update_iter_summary(self,n_iter,summary):
        with self.lock:
            self.we_h5file['summary'][n_iter-1] = summary

    def del_iter_summary(self, min_iter): #delete the iterations starting at min_iter      
        with self.lock:
            self.we_h5file['summary'].resize((min_iter - 1,))
                     
    def write_recycling_data(self, n_iter, rec_summary):
        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            rec_data_ds = iter_group.require_dataset('recycling', (len(rec_summary),), dtype=rec_summary_dtype)
            rec_data = numpy.zeros((len(rec_summary),), dtype=rec_summary_dtype)
            for itarget, target in enumerate(rec_summary):
                count, weight = target
                rec_data[itarget]['count'] = count
                rec_data[itarget]['weight'] = weight
            rec_data_ds[...] = rec_data
        
    def write_bin_data(self, n_iter, assignments, populations, n_trans, fluxes, rates):
        with self.lock:
            try:
                iter_group = self.get_iter_group(n_iter)
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
          * initpoint_type
        """
        
        segments = sorted(segments, key=attrgetter('seg_id'))
        
        if len(segments) == 1:
            self._update_segment(n_iter, segments[0])
            return
        
        seg_ids = [segment.seg_id for segment in segments]

        with self.lock:        
            iter_group = self.get_iter_group(n_iter)
            seg_index_entries = iter_group['seg_index'][seg_ids]
            pcoord_entries = iter_group['pcoord'][seg_ids]
            
            assert len(seg_index_entries) == len(pcoord_entries) == len(seg_ids)
                        
            row = numpy.empty((1,), seg_index_dtype)
            for (iseg, segment) in enumerate(segments):
                row[:] = seg_index_entries[iseg]
                row['status'] = segment.status
                row['endpoint_type'] = segment.endpoint_type or Segment.SEG_ENDPOINT_UNSET
                row['cputime'] = segment.cputime
                row['walltime'] = segment.walltime
                row['weight'] = segment.weight
                
                seg_index_entries[iseg] = row
                
                #pcoords[seg_id] = segment.pcoord
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
            
            
    def _update_segment(self, n_iter, segment):
        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            seg_index = self.get_seg_index(n_iter)
            row = seg_index[segment.seg_id]
            row['status'] = segment.status
            row['endpoint_type'] = segment.endpoint_type or Segment.SEG_ENDPOINT_UNSET
            row['cputime'] = segment.cputime
            row['walltime'] = segment.walltime
            row['weight'] = segment.weight
            seg_index[segment.seg_id] = row
            iter_group['pcoord'][segment.seg_id] = segment.pcoord
            
            if segment.data:
                aux_group = iter_group.require_group('auxdata')
                for (dsname, data) in segment.data.iteritems():
                    adata = numpy.asarray(data)
                    shape = (len(seg_index),)+adata.shape
                    nbytes = numpy.multiply.reduce(shape)
                    dset = aux_group.require_dataset(dsname,
                                                     shape=shape,
                                                     dtype=adata.dtype,
                                                     compression='gzip' if nbytes > self.aux_compression_threshold else None)
                    dset[segment.seg_id] = adata
            
                
    def get_segments(self, n_iter, load_auxdata=None):
        '''Return the segments from a given iteration.  This function is optimized for the 
        case of retrieving (nearly) all segments for a given iteration as quickly as possible, 
        and as such effectively loads all data for the given iteration into memory (which
        is what is currently required for running a WE iteration).'''
        
        if load_auxdata is None: load_auxdata = self.auto_load_auxdata
        
        with self.lock:
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
                                  initpoint_type = row['initpoint_type'],
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
        
        with self.lock:
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
                                  initpoint_type = row['initpoint_type'],
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
        
        with self.lock:
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
            # the above will return a scalar if only one is found, so convert
            # to a list if necessary
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
        
