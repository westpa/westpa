"""
HDF5 data manager for WEST.

Original HDF5 implementation: Joseph W. Kaus
Current implementation: Matthew C. Zwier

WEST exclusively uses the cross-platform, self-describing file format HDF5
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
            - auxdata/ -- auxiliary datasets (data stored on the 'data' field of Segment objects)

The file root object has an integer attribute 'west_file_format_version' which can be used to
determine how to access data even as the file format (i.e. organization of data within HDF5 file)
evolves.

Version history:
    Version 9
        - Basis states are now saved as iter_segid instead of just segid as a pointer label.
        - Initial states are also saved in the iteration 0 file, with a negative sign.
    Version 8
        - Added external links to trajectory files in iterations/iter_* groups, if the HDF5
          framework was used.
        - Added an iter group for the iteration 0 to store conformations of basis states.
    Version 7
        - Removed bin_assignments, bin_populations, and bin_rates from iteration group.
        - Added new_segments subgroup to iteration group
    Version 6
        - ???
    Version 5
        - moved iter_* groups into a top-level iterations/ group,
        - added in-HDF5 storage for basis states, target states, and generated states
"""

import logging
import pickle
import posixpath
import sys
import threading
import time
from operator import attrgetter
from os.path import relpath, dirname

import h5py
from h5py import h5s
import numpy as np

from . import h5io
from .segment import Segment
from .states import BasisState, TargetState, InitialState
from .we_driver import NewWeightEntry
from .propagators.executable import ExecutablePropagator

import westpa


log = logging.getLogger(__name__)

file_format_version = 9

makepath = ExecutablePropagator.makepath


class flushing_lock:
    def __init__(self, lock, fileobj):
        self.lock = lock
        self.fileobj = fileobj

    def __enter__(self):
        self.lock.acquire()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fileobj.flush()
        self.lock.release()


class expiring_flushing_lock:
    def __init__(self, lock, flush_method, nextsync):
        self.lock = lock
        self.flush_method = flush_method
        self.nextsync = nextsync

    def __enter__(self):
        self.lock.acquire()

    def __exit__(self, exc_type, exc_val, exc_tb):
        if time.time() > self.nextsync:
            self.flush_method()
        self.lock.release()


# Data types for use in the HDF5 file
seg_id_dtype = np.int64  # Up to 9 quintillion segments per iteration; signed so that initial states can be stored negative
n_iter_dtype = np.uint32  # Up to 4 billion iterations
weight_dtype = np.float64  # about 15 digits of precision in weights
utime_dtype = np.float64  # ("u" for Unix time) Up to ~10^300 cpu-seconds
vstr_dtype = h5py.special_dtype(vlen=str)
h5ref_dtype = h5py.special_dtype(ref=h5py.Reference)
binhash_dtype = np.dtype('|S64')

# seg_status_dtype    = h5py.special_dtype(enum=(np.uint8, Segment.statuses))
# seg_initpoint_dtype = h5py.special_dtype(enum=(np.uint8, Segment.initpoint_types))
# seg_endpoint_dtype  = h5py.special_dtype(enum=(np.uint8, Segment.endpoint_types))
# istate_type_dtype   = h5py.special_dtype(enum=(np.uint8, InitialState.istate_types))
# istate_status_dtype = h5py.special_dtype(enum=(np.uint8, InitialState.istate_statuses))

seg_status_dtype = np.uint8
seg_initpoint_dtype = np.uint8
seg_endpoint_dtype = np.uint8
istate_type_dtype = np.uint8
istate_status_dtype = np.uint8

summary_table_dtype = np.dtype(
    [
        ('n_particles', seg_id_dtype),  # Number of live trajectories in this iteration
        ('norm', weight_dtype),  # Norm of probability, to watch for errors or drift
        ('min_bin_prob', weight_dtype),  # Per-bin minimum probability
        ('max_bin_prob', weight_dtype),  # Per-bin maximum probability
        ('min_seg_prob', weight_dtype),  # Per-segment minimum probability
        ('max_seg_prob', weight_dtype),  # Per-segment maximum probability
        ('cputime', utime_dtype),  # Total CPU time for this iteration
        ('walltime', utime_dtype),  # Total wallclock time for this iteration
        ('binhash', binhash_dtype),
    ]
)


# The HDF5 file tracks two distinct, but related, histories:
#    (1) the evolution of the trajectory, which requires only an identifier
#        of where a segment's initial state comes from (the "history graph");
#        this is stored as the parent_id field of the seg index
#    (2) the flow of probability due to splits, merges, and recycling events,
#        which can be thought of as an adjacency list (the "weight graph")
# segment ID is implied by the row in the index table, and so is not stored
# initpoint_type remains implicitly stored as negative IDs (if parent_id < 0, then init_state_id = -(parent_id+1)
seg_index_dtype = np.dtype(
    [
        ('weight', weight_dtype),  # Statistical weight of this segment
        ('parent_id', seg_id_dtype),  # ID of parent (for trajectory history)
        ('wtg_n_parents', np.uint),  # number of parents this segment has in the weight transfer graph
        ('wtg_offset', np.uint),  # offset into the weight transfer graph dataset
        ('cputime', utime_dtype),  # CPU time used in propagating this segment
        ('walltime', utime_dtype),  # Wallclock time used in propagating this segment
        ('endpoint_type', seg_endpoint_dtype),  # Endpoint type (will continue, merged, or recycled)
        ('status', seg_status_dtype),  # Status of propagation of this segment
    ]
)

# Index to basis/initial states
ibstate_index_dtype = np.dtype([('iter_valid', np.uint), ('n_bstates', np.uint), ('group_ref', h5ref_dtype)])

# Basis state index type
bstate_dtype = np.dtype(
    [
        ('label', vstr_dtype),  # An optional descriptive label
        ('probability', weight_dtype),  # Probability that this state will be selected
        ('auxref', vstr_dtype),  # An optional auxiliar data reference
    ]
)

# Even when initial state generation is off and basis states are passed through directly, an initial state entry
# is created, as that allows precise tracing of the history of a given state in the most complex case of
# a new initial state for every new trajectory.
istate_dtype = np.dtype(
    [
        ('iter_created', np.uint),  # Iteration during which this state was generated (0 for at w_init)
        ('iter_used', np.uint),  # When this state was used to start a new trajectory
        ('basis_state_id', seg_id_dtype),  # Which basis state this state was generated from
        ('istate_type', istate_type_dtype),  # What type this initial state is (generated or basis)
        ('istate_status', istate_status_dtype),  # Whether this initial state is ready to go
        ('basis_auxref', vstr_dtype),
    ]
)

tstate_index_dtype = np.dtype(
    [('iter_valid', np.uint), ('n_states', np.uint), ('group_ref', h5ref_dtype)]  # Iteration when this state list is valid
)  # Reference to a group containing further data; this will be the
# null reference if there is no target state for that timeframe.
tstate_dtype = np.dtype([('label', vstr_dtype)])  # An optional descriptive label for this state

# Support for west.we_driver.NewWeightEntry
nw_source_dtype = np.uint8
nw_index_dtype = np.dtype(
    [
        ('source_type', nw_source_dtype),
        ('weight', weight_dtype),
        ('prev_seg_id', seg_id_dtype),
        ('target_state_id', seg_id_dtype),
        ('initial_state_id', seg_id_dtype),
    ]
)

# Storage of bin identities
binning_index_dtype = np.dtype([('hash', binhash_dtype), ('pickle_len', np.uint32)])


class WESTDataManager:
    """Data manager for assisiting the reading and writing of WEST data from/to HDF5 files."""

    # defaults for various options
    default_iter_prec = 8
    default_we_h5filename = 'west.h5'
    default_we_h5file_driver = None
    default_flush_period = 60

    # Compress any auxiliary dataset whose total size (across all segments) is more than 1MB
    default_aux_compression_threshold = 1048576

    # Bin data horizontal (second dimension) chunk size
    binning_hchunksize = 4096

    # Number of rows to retrieve during a table scan
    table_scan_chunksize = 1024

    def flushing_lock(self):
        return flushing_lock(self.lock, self.we_h5file)

    def expiring_flushing_lock(self):
        next_flush = self.last_flush + self.flush_period
        return expiring_flushing_lock(self.lock, self.flush_backing, next_flush)

    def process_config(self):
        config = self.rc.config

        for entry, type_ in [('iter_prec', int)]:
            config.require_type_if_present(['west', 'data', entry], type_)

        self.we_h5filename = config.get_path(['west', 'data', 'west_data_file'], default=self.default_we_h5filename)
        self.we_h5file_driver = config.get_choice(
            ['west', 'data', 'west_data_file_driver'],
            [None, 'sec2', 'family'],
            default=self.default_we_h5file_driver,
            value_transform=(lambda x: x.lower() if x else None),
        )
        self.iter_prec = config.get(['west', 'data', 'iter_prec'], self.default_iter_prec)
        self.aux_compression_threshold = config.get(
            ['west', 'data', 'aux_compression_threshold'], self.default_aux_compression_threshold
        )
        self.flush_period = config.get(['west', 'data', 'flush_period'], self.default_flush_period)
        self.iter_ref_h5_template = config.get(['west', 'data', 'data_refs', 'iteration'], None)
        self.store_h5 = self.iter_ref_h5_template is not None

        # Process dataset options
        dsopts_list = config.get(['west', 'data', 'datasets']) or []
        for dsopts in dsopts_list:
            dsopts = normalize_dataset_options(dsopts, path_prefix='auxdata' if dsopts['name'] != 'pcoord' else '')
            try:
                self.dataset_options[dsopts['name']].update(dsopts)
            except KeyError:
                self.dataset_options[dsopts['name']] = dsopts

        if 'pcoord' in self.dataset_options:
            if self.dataset_options['pcoord']['h5path'] != 'pcoord':
                raise ValueError('cannot override pcoord storage location')

    def __init__(self, rc=None):
        self.rc = rc or westpa.rc

        self.we_h5filename = self.default_we_h5filename
        self.we_h5file_driver = self.default_we_h5file_driver
        self.we_h5file_version = None
        self.h5_access_mode = 'r+'
        self.iter_prec = self.default_iter_prec
        self.aux_compression_threshold = self.default_aux_compression_threshold

        self.we_h5file = None

        self.lock = threading.RLock()
        self.flush_period = None
        self.last_flush = 0

        self._system = None
        self.iter_ref_h5_template = None
        self.store_h5 = False

        self.dataset_options = {}
        self.process_config()

    @property
    def system(self):
        if self._system is None:
            self._system = self.rc.get_system_driver()
        return self._system

    @system.setter
    def system(self, system):
        self._system = system

    @property
    def closed(self):
        return self.we_h5file is None

    def iter_group_name(self, n_iter, absolute=True):
        if absolute:
            return '/iterations/iter_{:0{prec}d}'.format(int(n_iter), prec=self.iter_prec)
        else:
            return 'iter_{:0{prec}d}'.format(int(n_iter), prec=self.iter_prec)

    def require_iter_group(self, n_iter):
        '''Get the group associated with n_iter, creating it if necessary.'''
        with self.lock:
            iter_group = self.we_h5file.require_group('/iterations/iter_{:0{prec}d}'.format(int(n_iter), prec=self.iter_prec))
            iter_group.attrs['n_iter'] = n_iter
        return iter_group

    def del_iter_group(self, n_iter):
        with self.lock:
            del self.we_h5file['/iterations/iter_{:0{prec}d}'.format(int(n_iter), prec=self.iter_prec)]

    def get_iter_group(self, n_iter):
        with self.lock:
            try:
                return self.we_h5file['/iterations/iter_{:0{prec}d}'.format(int(n_iter), prec=self.iter_prec)]
            except KeyError:
                return self.we_h5file['/iter_{:0{prec}d}'.format(int(n_iter), prec=self.iter_prec)]

    def get_seg_index(self, n_iter):
        with self.lock:
            seg_index = self.get_iter_group(n_iter)['seg_index']
            return seg_index

    @property
    def current_iteration(self):
        with self.lock:
            h5file_attrs = self.we_h5file['/'].attrs
            h5file_attr_keys = list(h5file_attrs.keys())

            if 'west_current_iteration' in h5file_attr_keys:
                return int(self.we_h5file['/'].attrs['west_current_iteration'])
            else:
                return int(self.we_h5file['/'].attrs['wemd_current_iteration'])

    @current_iteration.setter
    def current_iteration(self, n_iter):
        with self.lock:
            self.we_h5file['/'].attrs['west_current_iteration'] = n_iter

    def open_backing(self, mode=None):
        '''Open the (already-created) HDF5 file named in self.west_h5filename.'''
        mode = mode or self.h5_access_mode
        if not self.we_h5file:
            log.debug('attempting to open {} with mode {}'.format(self.we_h5filename, mode))
            self.we_h5file = h5io.WESTPAH5File(self.we_h5filename, mode, driver=self.we_h5file_driver)

            h5file_attrs = self.we_h5file['/'].attrs
            h5file_attr_keys = list(h5file_attrs.keys())

            if 'west_iter_prec' in h5file_attr_keys:
                self.iter_prec = int(h5file_attrs['west_iter_prec'])
            elif 'wemd_iter_prec' in h5file_attr_keys:
                self.iter_prec = int(h5file_attrs['wemd_iter_prec'])
            else:
                log.info('iteration precision not stored in HDF5; using {:d}'.format(self.iter_prec))

            if 'west_file_format_version' in h5file_attr_keys:
                self.we_h5file_version = h5file_attrs['west_file_format_version']
            elif 'wemd_file_format_version' in h5file_attr_keys:
                self.we_h5file_version = h5file_attrs['wemd_file_format_version']
            else:
                log.info('WEST HDF5 file format version not stored, assuming 0')
                self.we_h5file_version = 0

            log.debug('opened WEST HDF5 file version {:d}'.format(self.we_h5file_version))

    def prepare_backing(self):  # istates):
        '''Create new HDF5 file'''
        self.we_h5file = h5py.File(self.we_h5filename, 'w', driver=self.we_h5file_driver)

        with self.flushing_lock():
            self.we_h5file['/'].attrs['west_file_format_version'] = file_format_version
            self.we_h5file['/'].attrs['west_iter_prec'] = self.iter_prec
            self.we_h5file['/'].attrs['west_version'] = westpa.__version__
            self.current_iteration = 0
            self.we_h5file['/'].create_dataset('summary', shape=(1,), dtype=summary_table_dtype, maxshape=(None,))
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
                self.last_flush = time.time()

    def save_target_states(self, tstates, n_iter=None):
        '''Save the given target states in the HDF5 file; they will be used for the next iteration to
        be propagated.  A complete set is required, even if nominally appending to an existing set,
        which simplifies the mapping of IDs to the table.'''

        system = self.system

        n_iter = n_iter or self.current_iteration

        # Assemble all the important data before we start to modify the HDF5 file
        tstates = list(tstates)
        if tstates:
            state_table = np.empty((len(tstates),), dtype=tstate_dtype)
            state_pcoords = np.empty((len(tstates), system.pcoord_ndim), dtype=system.pcoord_dtype)
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
                master_index = master_group.create_dataset('index', shape=(1,), maxshape=(None,), dtype=tstate_index_dtype)
                n_sets = 1
            else:
                n_sets = len(master_index) + 1
                master_index.resize((n_sets,))

            set_id = n_sets - 1
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
            set_id = np.digitize([n_iter], master_index['iter_valid']) - 1
            group_ref = master_index[set_id]['group_ref']

            # Check if reference is Null
            if not bool(group_ref):
                return None

            # This extra [0] is to work around a bug in h5py
            try:
                group = self.we_h5file[group_ref]
            except (TypeError, AttributeError):
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

            if tstate_group is not None:
                tstate_index = tstate_group['index'][...]
                tstate_pcoords = tstate_group['pcoord'][...]

                tstates = [
                    TargetState(state_id=i, label=h5io.tostr(row['label']), pcoord=pcoord.copy())
                    for (i, (row, pcoord)) in enumerate(zip(tstate_index, tstate_pcoords))
                ]
            else:
                tstates = []

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
                master_index = master_group.create_dataset('index', dtype=ibstate_index_dtype, shape=(1,), maxshape=(None,))
                n_sets = 1
            else:
                n_sets = len(master_index) + 1
                master_index.resize((n_sets,))

            set_id = n_sets - 1
            master_index_row = master_index[set_id]
            master_index_row['iter_valid'] = n_iter
            master_index_row['n_bstates'] = len(basis_states)
            state_group = master_group.create_group(str(set_id))
            master_index_row['group_ref'] = state_group.ref

            if basis_states:
                system = self.system
                state_table = np.empty((len(basis_states),), dtype=bstate_dtype)
                state_pcoords = np.empty((len(basis_states), system.pcoord_ndim), dtype=system.pcoord_dtype)
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

    def create_ibstate_iter_h5file(self, basis_states):
        '''Create the per-iteration HDF5 file for the basis states (i.e., iteration 0).
        This special treatment is needed so that the analysis tools can access basis states
        more easily.'''

        if not self.store_h5:
            return

        segments = []
        for i, state in enumerate(basis_states):
            dummy_segment = Segment(
                n_iter=0,
                seg_id=state.state_id,
                parent_id=-(state.state_id + 1),
                weight=state.probability,
                wtg_parent_ids=None,
                pcoord=state.pcoord,
                status=Segment.SEG_STATUS_UNSET,
                data=state.data,
            )
            segments.append(dummy_segment)

        # # link the iteration file in west.h5
        self.prepare_iteration(0, segments)
        self.update_iter_h5file(0, segments)

    def update_iter_h5file(self, n_iter, segments):
        '''Write out the per-iteration HDF5 file with given segments and add an external link to it
        in the main HDF5 file (west.h5) if the link is not present.'''

        if not self.store_h5:
            return

        west_h5_file = makepath(self.we_h5filename)
        iter_ref_h5_file = makepath(self.iter_ref_h5_template, {'n_iter': n_iter})
        iter_ref_rel_path = relpath(iter_ref_h5_file, dirname(west_h5_file))

        with h5io.WESTIterationFile(iter_ref_h5_file, 'a') as outf:
            for segment in segments:
                outf.write_segment(segment, True)

        iter_group = self.get_iter_group(n_iter)

        if 'trajectories' not in iter_group:
            iter_group['trajectories'] = h5py.ExternalLink(iter_ref_rel_path, '/')

    def get_basis_states(self, n_iter=None):
        '''Return a list of BasisState objects representing the basis states that are in use for iteration n_iter.'''

        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)
            try:
                bstate_index = ibstate_group['bstate_index'][...]
            except KeyError:
                return []
            bstate_pcoords = ibstate_group['bstate_pcoord'][...]
            bstates = [
                BasisState(
                    state_id=i,
                    label=h5io.tostr(row['label']),
                    probability=row['probability'],
                    auxref=h5io.tostr(row['auxref']) or None,
                    pcoord=pcoord.copy(),
                )
                for (i, (row, pcoord)) in enumerate(zip(bstate_index, bstate_pcoords))
            ]

            bstate_total_prob = sum(bstate.probability for bstate in bstates)

            # This should run once in the second iteration, and only if start-states are specified,
            # but is necessary to re-normalize (i.e. normalize without start-state probabilities included)
            for i, bstate in enumerate(bstates):
                bstate.probability /= bstate_total_prob
                bstates[i] = bstate
            return bstates

    def create_initial_states(self, n_states, n_iter=None):
        '''Create storage for ``n_states`` initial states associated with iteration ``n_iter``, and
        return bare InitialState objects with only state_id set.'''

        system = self.system
        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)

            try:
                istate_index = ibstate_group['istate_index']
            except KeyError:
                istate_index = ibstate_group.create_dataset('istate_index', dtype=istate_dtype, shape=(n_states,), maxshape=(None,))
                istate_pcoords = ibstate_group.create_dataset(
                    'istate_pcoord',
                    dtype=system.pcoord_dtype,
                    shape=(n_states, system.pcoord_ndim),
                    maxshape=(None, system.pcoord_ndim),
                )
                len_index = len(istate_index)
                first_id = 0
            else:
                first_id = len(istate_index)
                len_index = len(istate_index) + n_states
                istate_index.resize((len_index,))
                istate_pcoords = ibstate_group['istate_pcoord']
                istate_pcoords.resize((len_index, system.pcoord_ndim))

        index_entries = istate_index[first_id:len_index]
        new_istates = []
        for irow, row in enumerate(index_entries):
            row['iter_created'] = n_iter
            row['istate_status'] = InitialState.ISTATE_STATUS_PENDING
            new_istates.append(
                InitialState(
                    state_id=first_id + irow,
                    basis_state_id=None,
                    iter_created=n_iter,
                    istate_status=InitialState.ISTATE_STATUS_PENDING,
                )
            )
        istate_index[first_id:len_index] = index_entries
        return new_istates

    def update_initial_states(self, initial_states, n_iter=None):
        '''Save the given initial states in the HDF5 file'''

        system = self.system
        initial_states = sorted(initial_states, key=attrgetter('state_id'))
        if not initial_states:
            return

        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)
            state_ids = [state.state_id for state in initial_states]
            index_entries = ibstate_group['istate_index'][state_ids]
            pcoord_vals = np.empty((len(initial_states), system.pcoord_ndim), dtype=system.pcoord_dtype)
            for i, initial_state in enumerate(initial_states):
                index_entries[i]['iter_created'] = initial_state.iter_created
                index_entries[i]['iter_used'] = initial_state.iter_used or InitialState.ISTATE_UNUSED
                index_entries[i]['basis_state_id'] = (
                    initial_state.basis_state_id if initial_state.basis_state_id is not None else -1
                )
                index_entries[i]['istate_type'] = initial_state.istate_type or InitialState.ISTATE_TYPE_UNSET
                index_entries[i]['istate_status'] = initial_state.istate_status or InitialState.ISTATE_STATUS_PENDING
                pcoord_vals[i] = initial_state.pcoord

                index_entries[i]['basis_auxref'] = initial_state.basis_auxref or ""

            ibstate_group['istate_index'][state_ids] = index_entries
            ibstate_group['istate_pcoord'][state_ids] = pcoord_vals

            if self.store_h5:
                segments = []
                for i, state in enumerate(initial_states):
                    dummy_segment = Segment(
                        n_iter=-state.iter_created,
                        seg_id=state.state_id,
                        parent_id=state.basis_state_id,
                        wtg_parent_ids=None,
                        pcoord=state.pcoord,
                        status=Segment.SEG_STATUS_PREPARED,
                        data=state.data,
                    )
                    segments.append(dummy_segment)
                self.update_iter_h5file(0, segments)

    def get_initial_states(self, n_iter=None):
        states = []
        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.find_ibstate_group(n_iter)
            try:
                istate_index = ibstate_group['istate_index'][...]
            except KeyError:
                return []
            istate_pcoords = ibstate_group['istate_pcoord'][...]

            for state_id, (state, pcoord) in enumerate(zip(istate_index, istate_pcoords)):
                states.append(
                    InitialState(
                        state_id=state_id,
                        basis_state_id=int(state['basis_state_id']),
                        iter_created=int(state['iter_created']),
                        iter_used=int(state['iter_used']),
                        istate_type=int(state['istate_type']),
                        basis_auxref=h5io.tostr(state['basis_auxref']),
                        pcoord=pcoord.copy(),
                    )
                )
            return states

    def get_segment_initial_states(self, segments, n_iter=None):
        '''Retrieve all initial states referenced by the given segments.'''

        with self.lock:
            n_iter = n_iter or self.current_iteration
            ibstate_group = self.get_iter_group(n_iter)['ibstates']

            istate_ids = {-int(segment.parent_id + 1) for segment in segments if segment.parent_id < 0}
            sorted_istate_ids = sorted(istate_ids)
            if not sorted_istate_ids:
                return []

            istate_rows = ibstate_group['istate_index'][sorted_istate_ids][...]
            istate_pcoords = ibstate_group['istate_pcoord'][sorted_istate_ids][...]
            istates = []

            for state_id, state, pcoord in zip(sorted_istate_ids, istate_rows, istate_pcoords):
                try:
                    b_auxref = h5io.tostr(state['basis_auxref'])
                except ValueError:
                    b_auxref = ''
                istate = InitialState(
                    state_id=state_id,
                    basis_state_id=int(state['basis_state_id']),
                    iter_created=int(state['iter_created']),
                    iter_used=int(state['iter_used']),
                    istate_type=int(state['istate_type']),
                    basis_auxref=b_auxref,
                    pcoord=pcoord.copy(),
                )
                istates.append(istate)
            return istates

    def get_unused_initial_states(self, n_states=None, n_iter=None):
        '''Retrieve any prepared but unused initial states applicable to the given iteration.
        Up to ``n_states`` states are returned; if ``n_states`` is None, then all unused states
        are returned.'''

        n_states = n_states or sys.maxsize
        ISTATE_UNUSED = InitialState.ISTATE_UNUSED
        ISTATE_STATUS_PREPARED = InitialState.ISTATE_STATUS_PREPARED
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
                istop = min(istart + chunksize, n_index_entries)
                istate_chunk = istate_index[istart:istop]
                pcoord_chunk = istate_pcoords[istart:istop]
                # state_ids = np.arange(istart,istop,dtype=np.uint)

                for ci in range(len(istate_chunk)):
                    row = istate_chunk[ci]
                    pcoord = pcoord_chunk[ci]
                    state_id = istart + ci
                    if row['iter_used'] == ISTATE_UNUSED and row['istate_status'] == ISTATE_STATUS_PREPARED:
                        istate = InitialState(
                            state_id=state_id,
                            basis_state_id=int(row['basis_state_id']),
                            iter_created=int(row['iter_created']),
                            iter_used=0,
                            istate_type=int(row['istate_type']),
                            pcoord=pcoord.copy(),
                            istate_status=ISTATE_STATUS_PREPARED,
                        )
                        states.append(istate)
                    del row, pcoord, state_id
                istart += chunksize
                del istate_chunk, pcoord_chunk  # , state_ids, unused, ids_of_unused
            log.debug('found {:d} unused states'.format(len(states)))
            return states[:n_states]

    def prepare_iteration(self, n_iter, segments):
        """Prepare for a new iteration by creating space to store the new iteration's data.
        The number of segments, their IDs, and their lineage must be determined and included
        in the set of segments passed in."""

        log.debug('preparing HDF5 group for iteration %d (%d segments)' % (n_iter, len(segments)))

        # Ensure we have a list for guaranteed ordering
        init = n_iter == 0
        segments = list(segments)
        n_particles = len(segments)
        system = self.system
        pcoord_ndim = system.pcoord_ndim
        pcoord_len = 2 if init else system.pcoord_len
        pcoord_dtype = system.pcoord_dtype

        with self.lock:
            if not init:
                # Create a table of summary information about each iteration
                summary_table = self.we_h5file['summary']
                if len(summary_table) < n_iter:
                    summary_table.resize((n_iter + 1,))

            iter_group = self.require_iter_group(n_iter)

            for linkname in ('seg_index', 'pcoord', 'wtgraph'):
                try:
                    del iter_group[linkname]
                except KeyError:
                    pass

            # everything indexed by [particle] goes in an index table
            seg_index_table_ds = iter_group.create_dataset('seg_index', shape=(n_particles,), dtype=seg_index_dtype)
            # unfortunately, h5py doesn't like in-place modification of individual fields; it expects
            # tuples. So, construct everything in a numpy array and then dump the whole thing into hdf5
            # In fact, this appears to be an h5py best practice (collect as much in ram as possible and then dump)
            seg_index_table = seg_index_table_ds[...]

            if not init:
                summary_row = np.zeros((1,), dtype=summary_table_dtype)
                summary_row['n_particles'] = n_particles
                summary_row['norm'] = np.add.reduce(list(map(attrgetter('weight'), segments)))
                summary_table[n_iter - 1] = summary_row

            # pcoord is indexed as [particle, time, dimension]
            pcoord_opts = self.dataset_options.get('pcoord', {'name': 'pcoord', 'h5path': 'pcoord', 'compression': False})
            shape = (n_particles, pcoord_len, pcoord_ndim)
            pcoord_ds = create_dataset_from_dsopts(iter_group, pcoord_opts, shape, pcoord_dtype)
            pcoord = np.empty((n_particles, pcoord_len, pcoord_ndim), pcoord_dtype)

            total_parents = 0
            for seg_id, segment in enumerate(segments):
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
                    if init:
                        if segment.pcoord.shape != pcoord.shape[2:]:
                            raise ValueError(
                                'basis state pcoord shape [%r] does not match expected shape [%r]'
                                % (segment.pcoord.shape, pcoord.shape[2:])
                            )
                        # Initial pcoord
                        pcoord[seg_id, 1, :] = segment.pcoord[:]
                    else:
                        if len(segment.pcoord) == 1:
                            # Initial pcoord
                            pcoord[seg_id, 0, :] = segment.pcoord[0, :]
                        elif segment.pcoord.shape != pcoord.shape[1:]:
                            raise ValueError(
                                'segment pcoord shape [%r] does not match expected shape [%r]'
                                % (segment.pcoord.shape, pcoord.shape[1:])
                            )
                        else:
                            pcoord[seg_id, ...] = segment.pcoord

            if total_parents > 0:
                wtgraph_ds = iter_group.create_dataset('wtgraph', (total_parents,), seg_id_dtype, compression='gzip', shuffle=True)
                parents = np.empty((total_parents,), seg_id_dtype)

                for seg_id, segment in enumerate(segments):
                    offset = seg_index_table[seg_id]['wtg_offset']
                    extent = seg_index_table[seg_id]['wtg_n_parents']
                    parent_list = list(segment.wtg_parent_ids)
                    parents[offset : offset + extent] = parent_list[:]

                    assert set(parents[offset : offset + extent]) == set(segment.wtg_parent_ids)

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

            tstate_group = self.find_tstate_group(n_iter)
            if tstate_group is not None:
                iter_group['tstates'] = tstate_group

    def get_iter_summary(self, n_iter=None):
        n_iter = n_iter or self.current_iteration
        with self.lock:
            return self.we_h5file['summary'][n_iter - 1]

    def update_iter_summary(self, summary, n_iter=None):
        n_iter = n_iter or self.current_iteration
        with self.lock:
            self.we_h5file['summary'][n_iter - 1] = summary

    def del_iter_summary(self, min_iter):  # delete the iterations starting at min_iter
        with self.lock:
            self.we_h5file['summary'].resize((min_iter - 1,))

    def update_segments(self, n_iter, segments):
        '''Update segment information in the HDF5 file; all prior information for each
        ``segment`` is overwritten, except for parent and weight transfer information.'''

        segments = sorted(segments, key=attrgetter('seg_id'))

        with self.lock:
            iter_group = self.get_iter_group(n_iter)

            pc_dsid = iter_group['pcoord'].id
            si_dsid = iter_group['seg_index'].id

            seg_ids = [segment.seg_id for segment in segments]
            n_segments = len(segments)
            n_total_segments = si_dsid.shape[0]
            system = self.system
            pcoord_ndim = system.pcoord_ndim
            pcoord_len = system.pcoord_len
            pcoord_dtype = system.pcoord_dtype

            seg_index_entries = np.empty((n_segments,), dtype=seg_index_dtype)
            pcoord_entries = np.empty((n_segments, pcoord_len, pcoord_ndim), dtype=pcoord_dtype)

            pc_msel = h5s.create_simple(pcoord_entries.shape, (h5s.UNLIMITED,) * pcoord_entries.ndim)
            pc_msel.select_all()
            si_msel = h5s.create_simple(seg_index_entries.shape, (h5s.UNLIMITED,))
            si_msel.select_all()
            pc_fsel = pc_dsid.get_space()
            si_fsel = si_dsid.get_space()

            for iseg in range(n_segments):
                seg_id = seg_ids[iseg]
                op = h5s.SELECT_OR if iseg != 0 else h5s.SELECT_SET
                si_fsel.select_hyperslab((seg_id,), (1,), op=op)
                pc_fsel.select_hyperslab((seg_id, 0, 0), (1, pcoord_len, pcoord_ndim), op=op)

            # read summary data so that we have valud parent and weight transfer information
            si_dsid.read(si_msel, si_fsel, seg_index_entries)

            for iseg, (segment, ientry) in enumerate(zip(segments, seg_index_entries)):
                ientry['status'] = segment.status
                ientry['endpoint_type'] = segment.endpoint_type or Segment.SEG_ENDPOINT_UNSET
                ientry['cputime'] = segment.cputime
                ientry['walltime'] = segment.walltime
                ientry['weight'] = segment.weight

                pcoord_entries[iseg] = segment.pcoord

            # write progress coordinates and index using low level HDF5 functions for efficiency
            si_dsid.write(si_msel, si_fsel, seg_index_entries)
            pc_dsid.write(pc_msel, pc_fsel, pcoord_entries)

            # Now, to deal with auxiliary data
            # If any segment has any auxiliary data, then the aux dataset must spring into
            # existence. Each is named according to the name in segment.data, and has shape
            # (n_total_segs, ...) where the ... is the shape of the data in segment.data (and may be empty
            # in the case of scalar data) and dtype is taken from the data type of the data entry
            # compression is on by default for datasets that will be more than 1MiB

            # a mapping of data set name to (per-segment shape, data type) tuples
            dsets = {}

            # First we scan for presence, shape, and data type of auxiliary data sets
            for segment in segments:
                if segment.data:
                    for dsname in segment.data:
                        if dsname.startswith('iterh5/'):
                            continue
                        data = np.asarray(segment.data[dsname], order='C')
                        segment.data[dsname] = data
                        dsets[dsname] = (data.shape, data.dtype)

            # Then we iterate over data sets and store data
            if dsets:
                for dsname, (shape, dtype) in dsets.items():
                    # dset = self._require_aux_dataset(iter_group, dsname, n_total_segments, shape, dtype)
                    try:
                        dsopts = self.dataset_options[dsname]
                    except KeyError:
                        dsopts = normalize_dataset_options({'name': dsname}, path_prefix='auxdata')

                    shape = (n_total_segments,) + shape
                    dset = require_dataset_from_dsopts(
                        iter_group, dsopts, shape, dtype, autocompress_threshold=self.aux_compression_threshold, n_iter=n_iter
                    )
                    if dset is None:
                        # storage is suppressed
                        continue
                    for segment in segments:
                        try:
                            auxdataset = segment.data[dsname]
                        except KeyError:
                            pass
                        else:
                            source_rank = len(auxdataset.shape)
                            source_sel = h5s.create_simple(auxdataset.shape, (h5s.UNLIMITED,) * source_rank)
                            source_sel.select_all()
                            dest_sel = dset.id.get_space()
                            dest_sel.select_hyperslab((segment.seg_id,) + (0,) * source_rank, (1,) + auxdataset.shape)
                            dset.id.write(source_sel, dest_sel, auxdataset)
                    if 'delram' in list(dsopts.keys()):
                        del dsets[dsname]

            self.update_iter_h5file(n_iter, segments)

    def get_segments(self, n_iter=None, seg_ids=None, load_pcoords=True):
        '''Return the given (or all) segments from a given iteration.

        If the optional parameter ``load_auxdata`` is true, then all auxiliary datasets
        available are loaded and mapped onto the ``data`` dictionary of each segment. If
        ``load_auxdata`` is None, then use the default ``self.auto_load_auxdata``, which can
        be set by the option ``load_auxdata`` in the ``[data]`` section of ``west.cfg``. This
        essentially requires as much RAM as there is per-iteration auxiliary data, so this
        behavior is not on by default.'''

        n_iter = n_iter or self.current_iteration
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
                seg_ids = list(range(len(seg_index_ds)))
                seg_index_entries = seg_index_ds[...]
                if load_pcoords:
                    pcoord_entries = iter_group['pcoord'][...]

            segments = []

            for iseg, (seg_id, row) in enumerate(zip(seg_ids, seg_index_entries)):
                segment = Segment(
                    seg_id=seg_id,
                    n_iter=n_iter,
                    status=int(row['status']),
                    endpoint_type=int(row['endpoint_type']),
                    walltime=float(row['walltime']),
                    cputime=float(row['cputime']),
                    weight=float(row['weight']),
                )

                if load_pcoords:
                    segment.pcoord = pcoord_entries[iseg]

                if file_version < 5:
                    wtg_n_parents = row['n_parents']
                    wtg_offset = row['parents_offset']
                    wtg_parent_ids = all_parent_ids[wtg_offset : wtg_offset + wtg_n_parents]
                    segment.parent_id = int(wtg_parent_ids[0])
                else:
                    wtg_n_parents = row['wtg_n_parents']
                    wtg_offset = row['wtg_offset']
                    wtg_parent_ids = all_parent_ids[wtg_offset : wtg_offset + wtg_n_parents]
                    segment.parent_id = int(row['parent_id'])
                segment.wtg_parent_ids = set(map(int, wtg_parent_ids))
                assert len(segment.wtg_parent_ids) == wtg_n_parents
                segments.append(segment)
            del all_parent_ids
            if load_pcoords:
                del pcoord_entries

            # If any other data sets are requested, load them as well
            for dsinfo in self.dataset_options.values():
                if dsinfo.get('load', False):
                    dsname = dsinfo['name']
                    try:
                        ds = iter_group[dsinfo['h5path']]
                    except KeyError:
                        ds = None

                    if ds is not None:
                        for segment in segments:
                            seg_id = segment.seg_id
                            segment.data[dsname] = ds[seg_id]

        return segments

    def prepare_segment_restarts(self, segments, basis_states=None, initial_states=None):
        '''Prepare the necessary folder and files given the data stored in parent per-iteration HDF5 file
        for propagating the simulation. ``basis_states`` and ``initial_states`` should be provided if the
        segments are newly created'''

        if not self.store_h5:
            return

        for segment in segments:
            if segment.parent_id < 0:
                if initial_states is None or basis_states is None:
                    raise ValueError('initial and basis states required for preparing the segments')
                initial_state = initial_states[segment.initial_state_id]
                # Check if it's a start state
                if initial_state.istate_type == InitialState.ISTATE_TYPE_START:
                    log.debug(
                        f'Skip reading start state file from per-iteration HDF5 file for initial state {segment.initial_state_id}'
                    )
                    continue
                else:
                    basis_state = basis_states[initial_state.basis_state_id]

                parent = Segment(n_iter=0, seg_id=basis_state.state_id)
            else:
                parent = Segment(n_iter=segment.n_iter - 1, seg_id=segment.parent_id)

            try:
                parent_iter_ref_h5_file = makepath(self.iter_ref_h5_template, {'n_iter': parent.n_iter})

                with h5io.WESTIterationFile(parent_iter_ref_h5_file, 'r') as outf:
                    outf.read_restart(parent)

                segment.data['iterh5/restart'] = parent.data['iterh5/restart']
            except Exception as e:
                print('could not prepare restart data for segment {}/{}: {}'.format(segment.n_iter, segment.seg_id, str(e)))

    def get_all_parent_ids(self, n_iter):
        file_version = self.we_h5file_version
        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']

            if file_version < 5:
                offsets = seg_index['parents_offset']
                all_parents = iter_group['parents'][...]
                return all_parents.take(offsets)
            else:
                return seg_index['parent_id']

    def get_parent_ids(self, n_iter, seg_ids=None):
        '''Return a sequence of the parent IDs of the given seg_ids.'''

        file_version = self.we_h5file_version

        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']

            if seg_ids is None:
                seg_ids = range(len(seg_index))

            if file_version < 5:
                offsets = seg_index['parents_offset']
                all_parents = iter_group['parents'][...]
                return [all_parents[offsets[seg_id]] for seg_id in seg_ids]
            else:
                all_parents = seg_index['parent_id']
                return [all_parents[seg_id] for seg_id in seg_ids]

    def get_weights(self, n_iter, seg_ids):
        '''Return the weights associated with the given seg_ids'''

        unique_ids = sorted(set(seg_ids))
        if not unique_ids:
            return []
        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            index_subset = iter_group['seg_index'][unique_ids]
            weight_map = dict(zip(unique_ids, index_subset['weight']))
            return [weight_map[seg_id] for seg_id in seg_ids]

    def get_child_ids(self, n_iter, seg_id):
        '''Return the seg_ids of segments who have the given segment as a parent.'''

        with self.lock:
            if n_iter == self.current_iteration:
                return []

            iter_group = self.get_iter_group(n_iter + 1)
            seg_index = iter_group['seg_index']
            seg_ids = np.arange(len(seg_index), dtype=seg_id_dtype)

            if self.we_h5file_version < 5:
                offsets = seg_index['parents_offset']
                all_parent_ids = iter_group['parents'][...]
                parent_ids = np.array([all_parent_ids[offset] for offset in offsets])
            else:
                parent_ids = seg_index['parent_id']

            return seg_ids[parent_ids == seg_id]

    def get_children(self, segment):
        '''Return all segments which have the given segment as a parent'''

        if segment.n_iter == self.current_iteration:
            return []

        # Examine the segment index from the following iteration to see who has this segment
        # as a parent.  We don't need to worry about the number of parents each segment
        # has, since each has at least one, and indexing on the offset into the parents array
        # gives the primary parent ID

        with self.lock:
            iter_group = self.get_iter_group(segment.n_iter + 1)
            seg_index = iter_group['seg_index'][...]

            # This is one of the slowest pieces of code I've ever written...
            # seg_index = iter_group['seg_index'][...]
            # seg_ids = [seg_id for (seg_id,row) in enumerate(seg_index)
            #           if all_parent_ids[row['parents_offset']] == segment.seg_id]
            # return self.get_segments_by_id(segment.n_iter+1, seg_ids)
            if self.we_h5file_version < 5:
                parents = iter_group['parents'][seg_index['parent_offsets']]
            else:
                parents = seg_index['parent_id']
            all_seg_ids = np.arange(seg_index.len(), dtype=np.uintp)
            seg_ids = all_seg_ids[parents == segment.seg_id]
            # the above will return a scalar if only one is found, so convert
            # to a list if necessary
            try:
                len(seg_ids)
            except TypeError:
                seg_ids = [seg_ids]

            return self.get_segments(segment.n_iter + 1, seg_ids)

    # The following are dictated by the SimManager interface
    def prepare_run(self):
        self.open_backing()

    def finalize_run(self):
        self.flush_backing()
        self.close_backing()

    def save_new_weight_data(self, n_iter, new_weights):
        '''Save a set of NewWeightEntry objects to HDF5. Note that this should
        be called for the iteration in which the weights appear in their
        new locations (e.g. for recycled walkers, the iteration following
        recycling).'''

        if not new_weights:
            return

        system = self.system

        index = np.empty(len(new_weights), dtype=nw_index_dtype)
        prev_init_pcoords = system.new_pcoord_array(len(new_weights))
        prev_final_pcoords = system.new_pcoord_array(len(new_weights))
        new_init_pcoords = system.new_pcoord_array(len(new_weights))

        for ientry, nwentry in enumerate(new_weights):
            row = index[ientry]
            row['source_type'] = nwentry.source_type
            row['weight'] = nwentry.weight
            row['prev_seg_id'] = nwentry.prev_seg_id
            # the following use -1 as a sentinel for a missing value
            row['target_state_id'] = nwentry.target_state_id if nwentry.target_state_id is not None else -1
            row['initial_state_id'] = nwentry.initial_state_id if nwentry.initial_state_id is not None else -1

            index[ientry] = row

            if nwentry.prev_init_pcoord is not None:
                prev_init_pcoords[ientry] = nwentry.prev_init_pcoord

            if nwentry.prev_final_pcoord is not None:
                prev_final_pcoords[ientry] = nwentry.prev_final_pcoord

            if nwentry.new_init_pcoord is not None:
                new_init_pcoords[ientry] = nwentry.new_init_pcoord

        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            try:
                del iter_group['new_weights']
            except KeyError:
                pass

            nwgroup = iter_group.create_group('new_weights')
            nwgroup['index'] = index
            nwgroup['prev_init_pcoord'] = prev_init_pcoords
            nwgroup['prev_final_pcoord'] = prev_final_pcoords
            nwgroup['new_init_pcoord'] = new_init_pcoords

    def get_new_weight_data(self, n_iter):
        with self.lock:
            iter_group = self.get_iter_group(n_iter)

            try:
                nwgroup = iter_group['new_weights']
            except KeyError:
                return []

            try:
                index = nwgroup['index'][...]
                prev_init_pcoords = nwgroup['prev_init_pcoord'][...]
                prev_final_pcoords = nwgroup['prev_final_pcoord'][...]
                new_init_pcoords = nwgroup['new_init_pcoord'][...]
            except (KeyError, ValueError):  # zero-length selections raise ValueError
                return []

        entries = []
        for i in range(len(index)):
            irow = index[i]

            prev_seg_id = irow['prev_seg_id']
            if prev_seg_id == -1:
                prev_seg_id = None

            initial_state_id = irow['initial_state_id']
            if initial_state_id == -1:
                initial_state_id = None

            target_state_id = irow['target_state_id']
            if target_state_id == -1:
                target_state_id = None

            entry = NewWeightEntry(
                source_type=irow['source_type'],
                weight=irow['weight'],
                prev_seg_id=prev_seg_id,
                prev_init_pcoord=prev_init_pcoords[i].copy(),
                prev_final_pcoord=prev_final_pcoords[i].copy(),
                new_init_pcoord=new_init_pcoords[i].copy(),
                target_state_id=target_state_id,
                initial_state_id=initial_state_id,
            )

            entries.append(entry)
        return entries

    def find_bin_mapper(self, hashval):
        '''Check to see if the given has value is in the binning table. Returns the index in the
        bin data tables if found, or raises KeyError if not.'''

        try:
            hashval = hashval.hexdigest()
        except AttributeError:
            pass

        with self.lock:
            # these will raise KeyError if the group doesn't exist, which also means
            # that bin data is not available, so no special treatment here
            try:
                binning_group = self.we_h5file['/bin_topologies']
                index = binning_group['index']
            except KeyError:
                raise KeyError('hash {} not found'.format(hashval))

            n_entries = len(index)
            if n_entries == 0:
                raise KeyError('hash {} not found'.format(hashval))

            chunksize = self.table_scan_chunksize
            for istart in range(0, n_entries, chunksize):
                chunk = index[istart : min(istart + chunksize, n_entries)]
                for i in range(len(chunk)):
                    if chunk[i]['hash'] == bytes(hashval, 'utf-8'):
                        return istart + i

            raise KeyError('hash {} not found'.format(hashval))

    def get_bin_mapper(self, hashval):
        '''Look up the given hash value in the binning table, unpickling and returning the corresponding
        bin mapper if available, or raising KeyError if not.'''

        # Convert to a hex digest if we need to
        try:
            hashval = hashval.hexdigest()
        except AttributeError:
            pass

        with self.lock:
            # these will raise KeyError if the group doesn't exist, which also means
            # that bin data is not available, so no special treatment here
            try:
                binning_group = self.we_h5file['/bin_topologies']
                index = binning_group['index']
                pkl = binning_group['pickles']
            except KeyError:
                raise KeyError('hash {} not found. Could not retrieve binning group'.format(hashval))

            n_entries = len(index)
            if n_entries == 0:
                raise KeyError('hash {} not found. No entries in index'.format(hashval))

            chunksize = self.table_scan_chunksize

            for istart in range(0, n_entries, chunksize):
                chunk = index[istart : min(istart + chunksize, n_entries)]
                for i in range(len(chunk)):
                    if chunk[i]['hash'] == bytes(hashval, 'utf-8'):
                        pkldat = bytes(pkl[istart + i, 0 : chunk[i]['pickle_len']].data)
                        mapper = pickle.loads(pkldat)
                        log.debug('loaded {!r} from {!r}'.format(mapper, binning_group))
                        log.debug('hash value {!r}'.format(hashval))
                        return mapper

            raise KeyError('hash {} not found'.format(hashval))

    def save_bin_mapper(self, hashval, pickle_data):
        '''Store the given mapper in the table of saved mappers. If the mapper cannot be stored,
        PickleError will be raised. Returns the index in the bin data tables where the mapper is stored.'''

        try:
            hashval = hashval.hexdigest()
        except AttributeError:
            pass
        pickle_data = bytes(pickle_data)

        # First, scan to see if the mapper already is in the HDF5 file
        try:
            return self.find_bin_mapper(hashval)
        except KeyError:
            pass

        # At this point, we have a valid pickle and know it's not stored
        with self.lock:
            binning_group = self.we_h5file.require_group('/bin_topologies')

            try:
                index = binning_group['index']
                pickle_ds = binning_group['pickles']
            except KeyError:
                index = binning_group.create_dataset('index', shape=(1,), maxshape=(None,), dtype=binning_index_dtype)
                pickle_ds = binning_group.create_dataset(
                    'pickles',
                    dtype=np.uint8,
                    shape=(1, len(pickle_data)),
                    maxshape=(None, None),
                    chunks=(1, 4096),
                    compression='gzip',
                    compression_opts=9,
                )
                n_entries = 1
            else:
                n_entries = len(index) + 1
                index.resize((n_entries,))
                new_hsize = max(pickle_ds.shape[1], len(pickle_data))
                pickle_ds.resize((n_entries, new_hsize))

            index_row = index[n_entries - 1]
            index_row['hash'] = hashval
            index_row['pickle_len'] = len(pickle_data)
            index[n_entries - 1] = index_row
            pickle_ds[n_entries - 1, : len(pickle_data)] = memoryview(pickle_data)
            return n_entries - 1

    def save_iter_binning(self, n_iter, hashval, pickled_mapper, target_counts):
        '''Save information about the binning used to generate segments for iteration n_iter.'''

        with self.lock:
            iter_group = self.get_iter_group(n_iter)

            try:
                del iter_group['bin_target_counts']
            except KeyError:
                pass

            iter_group['bin_target_counts'] = target_counts

            if hashval and pickled_mapper:
                self.save_bin_mapper(hashval, pickled_mapper)
                iter_group.attrs['binhash'] = hashval
            else:
                iter_group.attrs['binhash'] = ''


def normalize_dataset_options(dsopts, path_prefix='', n_iter=0):
    dsopts = dict(dsopts)

    ds_name = dsopts['name']
    if path_prefix:
        default_h5path = '{}/{}'.format(path_prefix, ds_name)
    else:
        default_h5path = ds_name

    dsopts.setdefault('h5path', default_h5path)
    dtype = dsopts.get('dtype')
    if dtype:
        if isinstance(dtype, str):
            dsopts['dtype'] = np.dtype(getattr(np, dtype))
        else:
            dsopts['dtype'] = np.dtype(dtype)

    dsopts['store'] = bool(dsopts['store']) if 'store' in dsopts else True
    dsopts['load'] = bool(dsopts['load']) if 'load' in dsopts else False

    return dsopts


def create_dataset_from_dsopts(group, dsopts, shape=None, dtype=None, data=None, autocompress_threshold=None, n_iter=None):
    # log.debug('create_dataset_from_dsopts(group={!r}, dsopts={!r}, shape={!r}, dtype={!r}, data={!r}, autocompress_threshold={!r})'
    #          .format(group,dsopts,shape,dtype,data,autocompress_threshold))
    if not dsopts.get('store', True):
        return None

    if 'file' in list(dsopts.keys()):
        import h5py

        #        dsopts['file'] = str(dsopts['file']).format(n_iter=n_iter)
        h5_auxfile = h5io.WESTPAH5File(dsopts['file'].format(n_iter=n_iter))
        h5group = group
        if not ("iter_" + str(n_iter).zfill(8)) in h5_auxfile:
            h5_auxfile.create_group("iter_" + str(n_iter).zfill(8))
        group = h5_auxfile[('/' + "iter_" + str(n_iter).zfill(8))]

    h5path = dsopts['h5path']
    containing_group_name = posixpath.dirname(h5path)
    h5_dsname = posixpath.basename(h5path)

    # ensure arguments are sane
    if not shape and data is None:
        raise ValueError('either shape or data must be provided')
    elif data is None and (shape and dtype is None):
        raise ValueError('both shape and dtype must be provided when data is not provided')
    elif shape and data is not None and not data.shape == shape:
        raise ValueError('explicit shape {!r} does not match data shape {!r}'.format(shape, data.shape))

    if data is not None:
        shape = data.shape
        if dtype is None:
            dtype = data.dtype
    # end argument sanity checks

    # figure out where to store this data
    if containing_group_name:
        containing_group = group.require_group(containing_group_name)
    else:
        containing_group = group

    # has user requested an explicit data type?
    # the extra np.dtype is an idempotent operation on true dtype
    # objects, but ensures that things like np.float32, which are
    # actually NOT dtype objects, become dtype objects
    h5_dtype = np.dtype(dsopts.get('dtype', dtype))

    compression = None
    scaleoffset = None
    shuffle = False

    # compress if 1) explicitly requested, or 2) dataset size exceeds threshold and
    # compression not explicitly prohibited
    compression_directive = dsopts.get('compression')
    if compression_directive is None:
        # No directive
        nbytes = np.multiply.reduce(shape) * h5_dtype.itemsize
        if autocompress_threshold and nbytes > autocompress_threshold:
            compression = 9
    elif compression_directive == 0:  # includes False
        # Compression prohibited
        compression = None
    else:  # compression explicitly requested
        compression = compression_directive

    # Is scale/offset requested?
    scaleoffset = dsopts.get('scaleoffset', None)
    if scaleoffset is not None:
        scaleoffset = int(scaleoffset)

    # We always shuffle if we compress (losslessly)
    if compression:
        shuffle = True
    else:
        shuffle = False

    need_chunks = any([compression, scaleoffset is not None, shuffle])

    # We use user-provided chunks if available
    chunks_directive = dsopts.get('chunks')
    if chunks_directive is None:
        chunks = None
    elif chunks_directive is True:
        chunks = calc_chunksize(shape, h5_dtype)
    elif chunks_directive is False:
        chunks = None
    else:
        chunks = tuple(chunks_directive[i] if chunks_directive[i] <= shape[i] else shape[i] for i in range(len(shape)))

    if not chunks and need_chunks:
        chunks = calc_chunksize(shape, h5_dtype)

    opts = {'shape': shape, 'dtype': h5_dtype, 'compression': compression, 'shuffle': shuffle, 'chunks': chunks}

    try:
        import h5py._hl.filters

        h5py._hl.filters._COMP_FILTERS['scaleoffset']
    except (ImportError, KeyError, AttributeError):
        # filter not available, or an unexpected version of h5py
        # use lossless compression instead
        opts['compression'] = True
    else:
        opts['scaleoffset'] = scaleoffset

    if log.isEnabledFor(logging.DEBUG):
        log.debug('requiring aux dataset {!r}, shape={!r}, opts={!r}'.format(h5_dsname, shape, opts))

    dset = containing_group.require_dataset(h5_dsname, **opts)

    if data is not None:
        dset[...] = data

    if 'file' in list(dsopts.keys()):
        import h5py

        if not dsopts['h5path'] in h5group:
            h5group[dsopts['h5path']] = h5py.ExternalLink(
                dsopts['file'].format(n_iter=n_iter), ("/" + "iter_" + str(n_iter).zfill(8) + "/" + dsopts['h5path'])
            )

    return dset


def require_dataset_from_dsopts(group, dsopts, shape=None, dtype=None, data=None, autocompress_threshold=None, n_iter=None):
    if not dsopts.get('store', True):
        return None
    try:
        return group[dsopts['h5path']]
    except KeyError:
        return create_dataset_from_dsopts(
            group, dsopts, shape=shape, dtype=dtype, data=data, autocompress_threshold=autocompress_threshold, n_iter=n_iter
        )


def calc_chunksize(shape, dtype, max_chunksize=262144):
    '''Calculate a chunk size for HDF5 data, anticipating that access will slice
    along lower dimensions sooner than higher dimensions.'''

    chunk_shape = list(shape)
    for idim in range(len(shape)):
        chunk_nbytes = np.multiply.reduce(chunk_shape) * dtype.itemsize
        while chunk_shape[idim] > 1 and chunk_nbytes > max_chunksize:
            chunk_shape[idim] >>= 1  # divide by 2
            chunk_nbytes = np.multiply.reduce(chunk_shape) * dtype.itemsize

        if chunk_nbytes <= max_chunksize:
            break

    chunk_shape = tuple(chunk_shape)
    log.debug(
        'selected chunk shape {} for data set of type {} shaped {} (chunk size = {} bytes)'.format(
            chunk_shape, dtype, shape, chunk_nbytes
        )
    )
    return chunk_shape
