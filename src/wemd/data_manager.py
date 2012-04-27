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
            - wtg_parents -- data used to reconstruct the split/merge history of trajectories
            - recycling -- flux and event count for recycled particles, on a per-target-state basis
            - aux_data/ -- auxiliary datasets (data stored on the 'data' field of Segment objects)

The file root object has an integer attribute 'wemd_file_format_version' which can be used to
determine how to access data even as the file format (i.e. organization of data within HDF5 file)
evolves. 

Version history:
    Version 5 
        - moved iter_* groups into a top-level iterations/ group,
        - added in-HDF5 storage for basis states, target states, and generated states
"""        
from __future__ import division; __metaclass__ = type
import warnings, sys
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
from wemd.states import BasisState, TargetState, InitialState

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
seg_id_dtype = numpy.int64  # Up to 9 quintillion segments per iteration; signed so that initial states can be stored negative
n_iter_dtype = numpy.uint32 # Up to 4 billion iterations
weight_dtype = numpy.float64  # about 15 digits of precision in weights
utime_dtype = numpy.float64  # ("u" for Unix time) Up to ~10^300 cpu-seconds 
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
istate_type_dtype = numpy.uint8
istate_status_dtype = numpy.uint8
    
summary_table_dtype = numpy.dtype( [ ('n_particles', seg_id_dtype),    # Number of live trajectories in this iteration
                                     ('norm', weight_dtype),          # Norm of probability, to watch for errors or drift
                                     ('min_bin_prob', weight_dtype),  # Per-bin minimum probability
                                     ('max_bin_prob', weight_dtype),  # Per-bin maximum probability
                                     ('min_seg_prob', weight_dtype),  # Per-segment minimum probability
                                     ('max_seg_prob', weight_dtype),  # Per-segment maximum probability
                                     ('cputime', utime_dtype),       # Total CPU time for this iteration
                                     ('walltime', utime_dtype),# Total wallclock time for this iteration
                                     ('binhash', '|S64'),
                                     ] )    


# The HDF5 file tracks two distinct, but related, histories: 
#    (1) the evolution of the trajectory, which requires only an identifier 
#        of where a segment's initial state comes from (the "history graph");
#        this is stored as the parent_id field of the seg index
#    (2) the flow of probability due to splits, merges, and recycling events,
#        which can be thought of as an adjacency list (the "weight graph")
# segment ID is implied by the row in the index table, and so is not stored
# initpoint_type remains implicitly stored as negative IDs (if parent_id < 0, then init_state_id = -(parent_id+1) 
seg_index_dtype = numpy.dtype( [ ('weight', weight_dtype),               # Statistical weight of this segment
                                 ('parent_id', seg_id_dtype),             # ID of parent (for trajectory history)
                                 ('wtg_n_parents', numpy.uint), # number of parents this segment has in the weight transfer graph
                                 ('wtg_offset', numpy.uint),    # offset into the weight transfer graph dataset
                                 ('cputime', utime_dtype),              # CPU time used in propagating this segment
                                 ('walltime', utime_dtype),             # Wallclock time used in propagating this segment
                                 ('endpoint_type', seg_endpoint_dtype),   # Endpoint type (will continue, merged, or recycled)
                                 ('status', seg_status_dtype),            # Status of propagation of this segment
                                 ] )

# Recycling summary, used for recording per-target flux
rec_summary_dtype = numpy.dtype( [ ('count', numpy.uint),        # Number of transitions into this state
                                   ('flux', weight_dtype) ] )   # Probability flux into this state

# Index to basis/initial states
ibstate_index_dtype = numpy.dtype([('iter_valid', numpy.uint),
                                   ('n_bstates', numpy.uint),
                                   ('group_ref', h5ref_dtype)])

# Basis state index type
bstate_dtype = numpy.dtype( [ ('label', vstr_dtype),            # An optional descriptive label
                              ('probability', weight_dtype),   # Probability that this state will be selected 
                              ('auxref', vstr_dtype),           # An optional auxiliar data reference 
                              ])

# Even when initial state generation is off and basis states are passed through directly, an initial state entry
# is created, as that allows precise tracing of the history of a given state in the most complex case of
# a new initial state for every new trajectory.
istate_dtype = numpy.dtype( [('iter_created', numpy.uint),      # Iteration during which this state was generated (0 for at w_init)
                             ('iter_used', numpy.uint),           # When this state was used to start a new trajectory 
                             ('basis_state_id', seg_id_dtype),    # Which basis state this state was generated from
                             ('istate_type', istate_type_dtype),  # What type this initial state is (generated or basis)
                             ('istate_status', istate_status_dtype), # Whether this initial state is ready to go
                             ]) 

tstate_index_dtype = numpy.dtype([('iter_valid', numpy.uint), # Iteration when this state list is valid
                                  ('n_states', numpy.uint),
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
        self.we_h5file_version = None
        self.h5_access_mode = 'r+'
        self.iter_prec = self.default_iter_prec
        self.aux_compression_threshold = self.default_aux_compression_threshold
        
        self.table_scan_chunksize = 256
        
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
            self.iter_prec      = wemd.rc.config.get_int('data.iter_prec', self.default_iter_prec)
            self.aux_compression_threshold = wemd.rc.config.get_int('data.aux_compression_threshold', 
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
            return '/iterations/iter_{:0{prec}d}'.format(long(n_iter), prec=self.iter_prec)
        else:
            return 'iter_{:0{prec}d}'.format(long(n_iter), prec=self.iter_prec)

    def require_iter_group(self, n_iter):
        '''Get the group associated with n_iter, creating it if necessary.'''
        with self.lock:
            iter_group = self.we_h5file.require_group('/iterations/iter_{:0{prec}d}'.format(long(n_iter), prec=self.iter_prec))
            iter_group.attrs['n_iter'] = n_iter
        return iter_group
            
    def del_iter_group(self, n_iter):
        with self.lock:
            del self.we_h5file['/iterations/iter_{:0{prec}d}'.format(long(n_iter), prec=self.iter_prec)]

    def get_iter_group(self, n_iter):
        with self.lock:
            try:
                return self.we_h5file['/iterations/iter_{:0{prec}d}'.format(long(n_iter), prec=self.iter_prec)]
            except KeyError:
                return self.we_h5file['/iter_{:0{prec}d}'.format(long(n_iter),prec=self.iter_prec)]
            
    def get_seg_index(self, n_iter):
        with self.lock:
            seg_index = self.get_iter_group(n_iter)['seg_index']
            return seg_index
        
    @property
    def current_iteration(self):
        with self.lock:
            return int(self.we_h5file['/'].attrs['wemd_current_iteration'])
    
    @current_iteration.setter
    def current_iteration(self, n_iter):
        with self.lock:
            self.we_h5file['/'].attrs['wemd_current_iteration'] = n_iter
        
    def open_backing(self, mode=None):
        '''Open the (already-created) HDF5 file named in self.wemd_h5filename.'''
        mode = mode or self.h5_access_mode
        if not self.we_h5file:
            log.debug('attempting to open {} with mode {}'.format(self.we_h5filename, mode))
            self.we_h5file = h5py.File(self.we_h5filename, mode, driver=self.we_h5file_driver)
            try:
                recorded_iter_prec = self.we_h5file['/'].attrs['wemd_iter_prec']
            except KeyError:
                log.info('iteration precision not stored in HDF5; using {:d}'.format(self.iter_prec))
            else:
                log.debug('iteration precision found: {:d}'.format(recorded_iter_prec))
                self.iter_prec = int(recorded_iter_prec)
                
            try:
                self.we_h5file_version = self.we_h5file['/'].attrs['wemd_file_format_version']
            except KeyError:
                log.info('WEMD HDF5 file format version not stored, assuming 0')
                self.we_h5file_version
            else:
                log.debug('opened WEMD HDF5 file version {:d}'.format(self.we_h5file_version))
                    
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

    def save_target_states(self, tstates, n_iter=None):
        '''Save the given target states in the HDF5 file; they will be used for the next iteration to
        be propagated.  A complete set is required, even if nominally appending to an existing set,
        which simplifies the mapping of IDs to the table.'''
        
        system = wemd.rc.get_system_driver()
        
        n_iter = n_iter or self.current_iteration

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
            master_group = self.we_h5file.require_group('tstates')
            
            try:
                master_index = master_group['index']
            except KeyError:
                master_index = master_group.create_dataset('index', shape=(1,), maxshape=(None,),
                                                            dtype=tstate_index_dtype)
                n_sets = 1
            else:
                n_sets = len(master_index) + 1
                master_index.resize((n_sets,))
            
            set_id = n_sets-1
            master_index_row = master_index[set_id]
            master_index_row['iter_valid'] = n_iter
            master_index_row['n_states'] = len(tstates)
            
            if tstates:
                state_group = master_group.create_group(str(set_id))
                master_index_row['group_ref'] = state_group.ref            
                state_group['index'] = state_table
                state_group['pcoord'] = state_pcoords
            else:
                master_index_row['group_ref'] = None 
                
            master_index[set_id] = master_index_row
            
    def _find_multi_iter_group(self, n_iter, master_group_name):
        with self.lock:
            master_group = self.we_h5file[master_group_name]
            master_index = master_group['index'][...]
            set_id = numpy.digitize([n_iter], master_index['iter_valid']) - 1
            # This extra [0] is to work around a bug in h5py
            group_ref = master_index[set_id]['group_ref']
            try:
                group = self.we_h5file[group_ref]
            except AttributeError:
                log.debug('working around h5py bug')
                group = self.we_h5file[group_ref[0]]
            else:
                log.debug('h5py fixed; remove alternate code path')
            log.debug('reference {!r} points to group {!r}'.format(group_ref, group))
            return group
            
    def find_tstate_group(self, n_iter):
        return self._find_multi_iter_group(n_iter, 'tstates')

    def find_ibstate_group(self, n_iter):
        return self._find_multi_iter_group(n_iter, 'ibstates')

    def get_target_states(self, n_iter):
        '''Return a list of Target objects representing the target (sink) states that are in use for iteration n_iter.
        Future iterations are assumed to continue from the most recent set of states.'''

        with self.lock:
            tstate_group = self.find_tstate_group(n_iter)         
            tstate_index = tstate_group['index']
            tstate_pcoords = tstate_group['pcoord']
            
            tstates = [TargetState(state_id=i, label=row['label'], pcoord=pcoord)
                        for (i, (row, pcoord))  in enumerate(izip(tstate_index, tstate_pcoords))]
            return tstates
          
    def create_ibstate_group(self, basis_states, n_iter=None):
        '''Create the group used to store basis states and initial states (whose definitions are always
        coupled).  This group is hard-linked into all iteration groups that use these basis and 
        initial states.'''
                
        with self.lock:
            n_iter = n_iter or self.current_iteration
            master_group = self.we_h5file.require_group('ibstates')
            
            try:
                master_index = master_group['index']
            except KeyError:
                master_index = master_group.create_dataset('index', dtype=ibstate_index_dtype,
                                                             shape=(1,), maxshape=(None,))
                n_sets = 1
            else:
                n_sets = len(master_index)+1
                master_index.resize((n_sets),)
            
            set_id = n_sets - 1
            master_index_row = master_index[set_id]
            master_index_row['iter_valid'] = n_iter
            master_index_row['n_bstates'] = len(basis_states)
            state_group = master_group.create_group(str(set_id))
            master_index_row['group_ref'] = state_group.ref
            
            
            if basis_states:
                system = wemd.rc.get_system_driver()            
                state_table = numpy.empty((len(basis_states),), dtype=bstate_dtype)
                state_pcoords = numpy.empty((len(basis_states),system.pcoord_ndim), dtype=system.pcoord_dtype)
                for i, state in enumerate(basis_states):
                    state.state_id = i
                    state_table[i]['label'] = state.label
                    state_table[i]['probability'] = state.probability
                    state_table[i]['auxref'] = state.auxref or ''
                    state_pcoords[i] = state.pcoord
                
                state_group['bstate_index'] = state_table
                state_group['bstate_pcoord'] = state_pcoords
            
            master_index[set_id] = master_index_row
            return state_group


    def get_basis_states(self, n_iter=None):
        '''Return a list of BasisState objects representing the basis states that are in use for iteration n_iter.'''
        
        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)
            bstate_index = ibstate_group['bstate_index']
            bstate_pcoords = ibstate_group['bstate_pcoord']
            bstates = [BasisState(state_id=i, label=row['label'], probability=row['probability'],
                                  auxref = row['auxref'] or None, pcoord=pcoord)
                       for (i, (row, pcoord))  in enumerate(izip(bstate_index, bstate_pcoords))]
            return bstates
            

    def create_initial_states(self, n_states, n_iter=None):
        '''Create storage for ``n_states`` initial states associated with iteration ``n_iter``, and
        return bare InitialState objects with only state_id set.'''
        
        system = wemd.rc.get_system_driver()        
        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)
            
            try:
                istate_index = ibstate_group['istate_index']
            except KeyError:
                istate_index = ibstate_group.create_dataset('istate_index', dtype=istate_dtype, 
                                                            shape=(n_states,), maxshape=(None,))
                istate_pcoords = ibstate_group.create_dataset('istate_pcoord', dtype=system.pcoord_dtype,
                                                              shape=(n_states,system.pcoord_ndim),
                                                              maxshape=(None,system.pcoord_ndim))
                len_index = len(istate_index)
                first_id = 0
            else:
                first_id = len(istate_index)
                len_index = len(istate_index) + n_states
                istate_index.resize((len_index,))
                istate_pcoords = ibstate_group['istate_pcoord']
                istate_pcoords.resize((len_index,system.pcoord_ndim))
                

        index_entries = istate_index[first_id:len_index]
        new_istates = []                
        for irow, row in enumerate(index_entries):
            row['iter_created'] = n_iter
            row['istate_status'] = InitialState.ISTATE_STATUS_PENDING
            new_istates.append(InitialState(state_id=first_id+irow, basis_state_id=None,
                                            iter_created=n_iter, istate_status=InitialState.ISTATE_STATUS_PENDING))
        istate_index[first_id:len_index] = index_entries
        return new_istates
            
    def update_initial_states(self, initial_states, n_iter = None):
        '''Save the given initial states in the HDF5 file'''
        
        system = wemd.rc.get_system_driver()
        initial_states = sorted(initial_states,key=attrgetter('state_id'))
        
        with self.lock:
            n_iter = n_iter or self.current_iteration     
            ibstate_group = self.find_ibstate_group(n_iter)
            state_ids = [state.state_id for state in initial_states]
            index_entries = ibstate_group['istate_index'][state_ids] 
            pcoord_vals = numpy.empty((len(initial_states), system.pcoord_ndim), dtype=system.pcoord_dtype)
            for i, initial_state in enumerate(initial_states):
                index_entries[i]['iter_created'] = initial_state.iter_created
                index_entries[i]['iter_used'] = initial_state.iter_used or InitialState.ISTATE_UNUSED
                index_entries[i]['basis_state_id'] = initial_state.basis_state_id
                index_entries[i]['istate_type'] = initial_state.istate_type or InitialState.ISTATE_TYPE_UNSET
                index_entries[i]['istate_status'] = initial_state.istate_status or InitialState.ISTATE_STATUS_PENDING
                pcoord_vals[i] = initial_state.pcoord
            
            ibstate_group['istate_index'][state_ids] = index_entries
            ibstate_group['istate_pcoord'][state_ids] = pcoord_vals
            
    def get_initial_states(self, n_iter=None):
        states = []
        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)
            istate_pcoords = ibstate_group['pcoord']
            istate_index = ibstate_group['istate_index']
            for state_id, (state, pcoord) in enumerate(izip(istate_index, istate_pcoords)):
                states.append(InitialState(state_id=state_id, basis_state_id=state['basis_state_id'],
                                           iter_created=state['iter_created'], iter_used=state['iter_used'],
                                           istate_type=state['istate_type'], pcoord=pcoord))
            return states
                

    def get_segment_initial_states(self, segments, n_iter=None):
        '''Retrieve all initial states referenced by the given segments.'''
        
        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.get_iter_group(n_iter)['ibstates']
            
            istate_ids = {-(segment.parent_id+1) for segment in segments if segment.parent_id < 0}
            sorted_istate_ids = sorted(istate_ids)
            if not sorted_istate_ids:
                return []
            
            istate_rows = ibstate_group['istate_index'][sorted_istate_ids]
            istate_pcoords = ibstate_group['istate_pcoord'][sorted_istate_ids]
            
            return [InitialState(state_id=state_id, basis_state_id=state['basis_state_id'],
                                 iter_created=state['iter_created'], iter_used=state['iter_used'],
                                 istate_type=state['istate_type'], pcoord=pcoord)
                    for state_id, state, pcoord in izip(sorted_istate_ids, istate_rows, istate_pcoords)] 
            
    def get_unused_initial_states(self, n_states = None, n_iter = None):
        '''Retrieve any prepared but unused initial states applicable to the given iteration.
        Up to ``n_states`` states are returned; if ``n_states`` is None, then all unused states
        are returned.'''
        n_states = n_states or sys.maxint
        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)
            istate_index = ibstate_group['istate_index']
            istate_pcoords = ibstate_group['istate_pcoord']
            n_index_entries = istate_index.len()
            chunksize = self.table_scan_chunksize
                        
            states = []
            istart = 0
            while istart < n_index_entries and len(states) < n_states:
                istop = min(istart+chunksize, n_index_entries)
                istate_chunk = istate_index[istart:istop]
                state_ids = numpy.arange(istart,istop,dtype=numpy.uint)
                
                unused = (  (istate_chunk['iter_used'] == InitialState.ISTATE_UNUSED) 
                          & (istate_chunk['istate_status'] == InitialState.ISTATE_STATUS_PREPARED ) )
                # TODO: add a check for invalid statuses and purge -- w_lint?
                ids_of_unused = list(state_ids[unused])
                if len(ids_of_unused):
                    pcoords = istate_pcoords[ids_of_unused]
                    states.extend(InitialState(state_id=state_id, basis_state_id=state['basis_state_id'],
                                               iter_created=state['iter_created'], iter_used=0, istate_type=state['istate_type'],
                                               pcoord=pcoord,
                                               istate_status=state['istate_status'])
                                  for state_id,state,pcoord in izip(ids_of_unused,istate_chunk[unused],pcoords))
                istart += chunksize
            log.debug('found {:d} unused states'.format(len(states)))
            return states[:n_states]
                
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
            
            iter_group = self.require_iter_group(n_iter)
            
            for linkname in ('seg_index', 'pcoord', 'wtgraph'):
                try:
                    del iter_group[linkname]
                except KeyError:
                    pass
            
            # Redundant, but this saves us having to parse the name to find the iteration number
            # if we iterate over all the groups
            iter_group.attrs['n_iter'] = n_iter
            
            # everything indexed by [particle] goes in an index table
            seg_index_table_ds = iter_group.create_dataset('seg_index', shape=(n_particles,),
                                                           dtype=seg_index_dtype)
            # unfortunately, h5py doesn't like in-place modification of individual fields; it expects
            # tuples. So, construct everything in a numpy array and then dump the whole thing into hdf5
            # In fact, this appears to be an h5py best practice (collect as much in ram as possible and then dump)
            seg_index_table = seg_index_table_ds[...]
                    
            summary_row = numpy.zeros((1,), dtype=summary_table_dtype)
            summary_row['n_particles'] = n_particles
            summary_row['norm'] = numpy.add.reduce(map(self._attrgetters['weight'], segments))
            summary_table[n_iter-1] = summary_row
            
            # pcoord is indexed as [particle, time, dimension]
            pcoord_ds = iter_group.create_dataset('pcoord', 
                                                  shape=(n_particles, pcoord_len, pcoord_ndim), 
                                                  dtype=pcoord_dtype)
            pcoord = numpy.empty((n_particles, pcoord_len, pcoord_ndim), pcoord_dtype)
                                    
            
            total_parents = 0
            for (seg_id, segment) in enumerate(segments):
                if segment.seg_id is not None:
                    assert segment.seg_id == seg_id
                else:
                    segment.seg_id = seg_id
                # Parent must be set, though what it means depends on initpoint_type
                assert segment.parent_id is not None
                segment.seg_id = seg_id
                seg_index_table[seg_id]['status'] = segment.status
                seg_index_table[seg_id]['weight'] = segment.weight
                seg_index_table[seg_id]['parent_id'] = segment.parent_id                
                seg_index_table[seg_id]['wtg_n_parents'] = len(segment.wtg_parent_ids)
                seg_index_table[seg_id]['wtg_offset'] = total_parents
                total_parents += len(segment.wtg_parent_ids)
    
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
                    
                        
            if total_parents > 0:
                wtgraph_ds = iter_group.create_dataset('wtgraph', (total_parents,), seg_id_dtype,
                                                       compression='gzip', shuffle=True)
                parents = numpy.empty((total_parents,), seg_id_dtype)
            
                for (seg_id, segment) in enumerate(segments):
                    offset = seg_index_table[seg_id]['wtg_offset']
                    extent = seg_index_table[seg_id]['wtg_n_parents']
                    parent_list = list(segment.wtg_parent_ids)
                    parents[offset:offset+extent] = parent_list[:]
                    
                    assert set(parents[offset:offset+extent]) == set(segment.wtg_parent_ids)
                
                wtgraph_ds[:] = parents

            # Create convenient hard links
            self.update_iter_group_links(n_iter)
            
    
            # Since we accumulated many of these changes in RAM (and not directly in HDF5), propagate
            # the changes out to HDF5
            seg_index_table_ds[:] = seg_index_table
            pcoord_ds[...] = pcoord

    def update_iter_group_links(self, n_iter):
        '''Update the per-iteration hard links pointing to the tables of target and initial/basis states for the
        given iteration.  These links are not used by this class, but are remarkably convenient for third-party
        analysis tools and hdfview.'''
        
        with self.lock:
            iter_group = self.require_iter_group(n_iter)
            
            for linkname in ('ibstates', 'tstates'):
                try:
                    del iter_group[linkname]
                except KeyError:
                    pass
            
            iter_group['ibstates'] = self.find_ibstate_group(n_iter)
            iter_group['tstates'] = self.find_tstate_group(n_iter)
            
    def get_iter_summary(self,n_iter=None):
        n_iter = n_iter or self.current_iteration
        with self.lock:
            return self.we_h5file['summary'][n_iter-1]
        
    def update_iter_summary(self,summary,n_iter=None):
        n_iter = n_iter or self.current_iteration
        with self.lock:
            self.we_h5file['summary'][n_iter-1] = summary

    def del_iter_summary(self, min_iter): #delete the iterations starting at min_iter      
        with self.lock:
            self.we_h5file['summary'].resize((min_iter - 1,))
                     
    def save_recycling_data(self, rec_summary, n_iter=None):
        with self.lock:
            n_iter = n_iter or self.current_iteration
            iter_group = self.get_iter_group(n_iter)
            rec_data_ds = iter_group.require_dataset('recycling', (len(rec_summary),), dtype=rec_summary_dtype)
            rec_data = rec_data_ds[...]
            for itarget, target in enumerate(rec_summary):
                count, flux = target
                rec_data[itarget]['count'] = count
                rec_data[itarget]['flux'] = flux
            rec_data_ds[...] = rec_data
        
    def save_bin_data(self, assignments, populations, n_trans, fluxes, rates, n_iter=None):
        with self.lock:
            n_iter = n_iter or self.current_iteration
            iter_group = self.get_iter_group(n_iter)
            for dataset in ('bin_assignments', 'bin_populations', 'bin_ntrans', 'bin_fluxes', 'bin_rates'):
                try:
                    del iter_group[dataset]
                except KeyError:
                    pass

            iter_group.create_dataset('bin_assignments', data=assignments, dtype=numpy.min_scalar_type(assignments.max()))                            
            iter_group.create_dataset('bin_populations', data=populations, dtype=weight_dtype)
            iter_group.create_dataset('bin_ntrans', data=n_trans, dtype=numpy.min_scalar_type(n_trans.max()))
            iter_group.create_dataset('bin_fluxes', data=fluxes, dtype=weight_dtype)
            iter_group.create_dataset('bin_rates', data=rates, dtype=weight_dtype)
        
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
                        
            for (iseg, (segment, ientry)) in enumerate(izip(segments,seg_index_entries)):
                ientry['status'] = segment.status
                ientry['endpoint_type'] = segment.endpoint_type or Segment.SEG_ENDPOINT_UNSET
                ientry['cputime'] = segment.cputime
                ientry['walltime'] = segment.walltime
                ientry['weight'] = segment.weight
                
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
    
    def get_segments(self, n_iter=None, seg_ids=None, load_pcoords = True, load_auxdata=None):
        '''Return the given (or all) segments from a given iteration.  
        
        If the optional parameter ``load_auxdata`` is true, then all auxiliary datasets
        available are loaded and mapped onto the ``data`` dictionary of each segment. If 
        ``load_auxdata`` is None, then use the default ``self.auto_load_auxdata``, which can
        be set by the option ``load_auxdata`` in the ``[data]`` section of ``wemd.cfg``. This
        essentially requires as much RAM as there is per-iteration auxiliary data, so this
        behavior is not on by default.'''
        

        n_iter = n_iter or self.current_iteration
        if load_auxdata is None: load_auxdata = self.auto_load_auxdata
        file_version = self.we_h5file_version
        
        with self.lock:            
            iter_group = self.get_iter_group(n_iter)
            seg_index_ds = iter_group['seg_index']
            
            if file_version < 5:
                all_parent_ids = iter_group['parents'][...]
            else:
                all_parent_ids = iter_group['wtgraph'][...]
            
            if seg_ids is not None:
                seg_ids = list(sorted(seg_ids))
                seg_index_entries = seg_index_ds[seg_ids]
                if load_pcoords:
                    pcoord_entries = iter_group['pcoord'][seg_ids]
            else:
                seg_ids = range(len(seg_index_ds))
                seg_index_entries = seg_index_ds[...]
                if load_pcoords:
                    pcoord_entries = iter_group['pcoord'][...]
                            
            segments = []
                
            for iseg, (seg_id, row) in enumerate(izip(seg_ids, seg_index_entries)):
                segment = Segment(seg_id = seg_id,
                                  n_iter = n_iter,
                                  status = int(row['status']),
                                  endpoint_type = int(row['endpoint_type']),
                                  walltime = float(row['walltime']),
                                  cputime = float(row['cputime']),
                                  weight = float(row['weight']))
                
                if load_pcoords:
                    segment.pcoord = pcoord_entries[iseg]
                    
                if file_version < 5:
                    wtg_n_parents = row['n_parents'] 
                    wtg_offset = row['parents_offset']
                    wtg_parent_ids = all_parent_ids[wtg_offset:wtg_offset+wtg_n_parents]
                    segment.parent_id = long(wtg_parent_ids[0])
                else:
                    wtg_n_parents = row['wtg_n_parents']
                    wtg_offset = row['wtg_offset']  
                    wtg_parent_ids = all_parent_ids[wtg_offset:wtg_offset+wtg_n_parents]
                    segment.parent_id = long(row['parent_id'])
                segment.wtg_parent_ids = set(imap(long,wtg_parent_ids))
                assert len(segment.wtg_parent_ids) == wtg_n_parents
                segments.append(segment)
            del all_parent_ids
            if load_pcoords:
                del pcoord_entries
        
            # If any auxiliary data sets are available, load them as well
            if load_auxdata and 'auxdata' in iter_group:
                for (dsname, ds) in iter_group['auxdata'].iteritems():
                    for (seg_id, segment) in enumerate(segments):
                        segment.data[dsname] = ds[seg_id]
        
        return segments
        
    def get_segments_by_id(self, n_iter, seg_ids, load_auxdata=None):
        warn_deprecated_usage('get_segments_by_id is deprecated; use get_segments(seg_ids=...) instead')
        return self.get_segments(n_iter, seg_ids, load_auxdata=load_auxdata)
    
    def get_parent_ids(self, n_iter, seg_ids=None):
        '''Return a sequence of the parent IDs of the given seg_ids.'''
        
        file_version = self.we_h5file_version
        
        
        with self.lock:
            iter_group = self.get_iter_group(n_iter)

            if seg_ids is not None:
                unique_ids = sorted(set(seg_ids))
                if not unique_ids:
                    return []
            else:
                seg_ids = unique_ids = range(iter_group['seg_index'].shape[0])            
            
            index_subset = iter_group['seg_index'][unique_ids]
            
            if file_version < 5:
                offsets = list(index_subset['parents_offset'])
                parent_map = dict(izip(unique_ids, iter_group['parents'][offsets]))                
            else:
                parent_map = dict(izip(unique_ids, index_subset['parent_id']))
                
        return [parent_map[seg_id] for seg_id in seg_ids]
    
    def get_weights(self, n_iter, seg_ids):
        '''Return the weights associated with the given seg_ids'''
        
        unique_ids = sorted(set(seg_ids))
        if not unique_ids:
            return []
        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            index_subset = iter_group['seg_index'][unique_ids]
            weight_map = dict(izip(unique_ids, index_subset['weight']))
            return [weight_map[seg_id] for seg_id in seg_ids]
                    
    def get_children(self, segment):
        '''Return all segments which have the given segment as a parent'''

        if segment.n_iter == self.current_iteration: return []
        
        # Examine the segment index from the following iteration to see who has this segment
        # as a parent.  We don't need to worry about the number of parents each segment
        # has, since each has at least one, and indexing on the offset into the parents array 
        # gives the primary parent ID
        
        with self.lock:
            iter_group = self.get_iter_group(segment.n_iter+1)
            seg_index = iter_group['seg_index'][...]
    
            # This is one of the slowest pieces of code I've ever written...
            #seg_index = iter_group['seg_index'][...]
            #seg_ids = [seg_id for (seg_id,row) in enumerate(seg_index) 
            #           if all_parent_ids[row['parents_offset']] == segment.seg_id]
            #return self.get_segments_by_id(segment.n_iter+1, seg_ids)
            if self.we_h5file_version < 5:
                parents = iter_group['parents'][seg_index['parent_offsets']] 
            else:
                parents = seg_index['parent_id']
            all_seg_ids = numpy.arange(seg_index.len(), dtype=numpy.uintp)
            seg_ids = all_seg_ids[parents == segment.seg_id]
            # the above will return a scalar if only one is found, so convert
            # to a list if necessary
            try:
                len(seg_ids)
            except TypeError:
                seg_ids = [seg_ids]
            
            return self.get_segments(segment.n_iter+1, seg_ids)

    # The following are dictated by the SimManager interface
    def prepare_run(self):
        self.open_backing()
                
    def finalize_run(self):
        self.flush_backing()
        self.close_backing()
        
