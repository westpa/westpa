import logging
from operator import attrgetter

import h5py
import numpy as np
from h5py import h5s

from westpa.core.data_manager import (
    WESTDataManager,
    seg_index_dtype,
    seg_id_dtype,
    summary_table_dtype,
    normalize_dataset_options,
    require_dataset_from_dsopts,
)
from westpa.core.state import State

logger = logging.getLogger(__name__)


def state_to_numpy(state):
    if state.coord is not None:
        dtype = np.dtype([('coord', state.coord.dtype, state.coord.shape)])
        return np.array(state.coord, dtype=dtype)
    else:
        dtype = np.dtype([('file', h5py.special_dtype(vlen=str))])
        return np.array(state.file, dtype=dtype)


def state_from_numpy(array):
    if 'coord' in array.dtype.names:
        return State(coord=array['coord'])
    else:
        return State(file=array['file'].decode('utf-8'))


class DataManager(WESTDataManager):

    def __init__(self, h5filename):
        super().__init__()
        self.we_h5filename = h5filename

    def prepare_iteration(self, n_iter, segments):
        logger.debug('preparing HDF5 group for iteration %d (%d segments)' % (n_iter, len(segments)))

        # Ensure we have a list for guaranteed ordering
        init = n_iter == 0
        segments = list(segments)
        n_particles = len(segments)

        with self.lock:
            if not init:
                # Create a table of summary information about each iteration
                summary_table = self.we_h5file['summary']
                if len(summary_table) < n_iter:
                    summary_table.resize((n_iter + 1,))

            iter_group = self.require_iter_group(n_iter)

            for dsname in ('seg_index', 'wtgraph', 'initial_states', 'final_states', 'pcoord'):
                try:
                    del iter_group[dsname]
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

            # Since we accumulated many of these changes in RAM (and not directly in HDF5), propagate
            # the changes out to HDF5
            seg_index_table_ds[:] = seg_index_table

            prepared_segments = [s for s in segments if s.initial_state is not None]
            if prepared_segments:
                self.write_initial_states(n_iter, prepared_segments)

    def finalize_iteration(self, n_iter, segments):
        self.update_seg_index(n_iter, segments)
        self.write_auxdata(n_iter, segments)

    def update_seg_index(self, n_iter, segments):
        """Update the ``seg_index`` dataset for a given iteration. All prior
        information for each segment is overwritten, except for parent and
        weight transfer information.

        Parameters
        ----------
        n_iter : int
            Iteration number.
        segments : iterable of Segment
            Set of segments belonging to iteration `n_iter`.

        """
        segments = sorted(segments, key=attrgetter('seg_id'))
        n_segments = len(segments)

        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            dsid = iter_group['seg_index'].id

            entries = np.empty((n_segments,), dtype=seg_index_dtype)

            msel = h5s.create_simple(entries.shape, (h5s.UNLIMITED,))
            msel.select_all()
            fsel = dsid.get_space()
            for i, segment in enumerate(segments):
                op = h5s.SELECT_OR if i != 0 else h5s.SELECT_SET
                fsel.select_hyperslab((segment.seg_id,), (1,), op=op)

            dsid.read(msel, fsel, entries)

            for segment, entry in zip(segments, entries):
                entry['status'] = segment.status
                entry['endpoint_type'] = segment.endpoint_type
                entry['cputime'] = segment.cputime
                entry['walltime'] = segment.walltime
                entry['weight'] = segment.weight

            dsid.write(msel, fsel, entries)

    def write_auxdata(self, n_iter, segments):
        # Now, to deal with auxiliary data
        # If any segment has any auxiliary data, then the aux dataset must spring into
        # existence. Each is named according to the name in segment.data, and has shape
        # (n_particles, ...) where the ... is the shape of the data in segment.data (and may be empty
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
        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            n_particles = iter_group['seg_index'].shape[0]

            if dsets:
                for dsname, (shape, dtype) in dsets.items():
                    try:
                        dsopts = self.dataset_options[dsname]
                    except KeyError:
                        dsopts = normalize_dataset_options({'name': dsname}, path_prefix='auxdata')

                    shape = (n_particles,) + shape
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

    def write_initial_states(self, n_iter, segments):
        """Write :attr:`~westpa.Segment.initial_state` data for a set of segments.

        Parametes
        ---------
        n_iter : int
            Iteration number.
        segments : list of Segment
            Set of segments belonging to iteration `n_iter`.

        """
        segments = sorted(segments, key=attrgetter('seg_id'))
        arrays = [state_to_numpy(s.initial_state) for s in segments]
        dtype = arrays[0].dtype  # infer dtype from first segment
        entries = np.fromiter(arrays, dtype=dtype)

        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            n_total_segments = iter_group['seg_index'].shape[0]

            ds = iter_group.require_dataset('initial_states', (n_total_segments,), dtype=dtype)
            dsid = ds.id

            msel = h5s.create_simple(entries.shape, (h5s.UNLIMITED,))
            msel.select_all()
            fsel = dsid.get_space()
            for i, segment in enumerate(segments):
                op = h5s.SELECT_OR if i != 0 else h5s.SELECT_SET
                fsel.select_hyperslab((segment.seg_id,), (1,), op=op)

            dsid.write(msel, fsel, entries)

    def write_final_states(self, n_iter, segments):
        """Write :attr:`~westpa.Segment.final_state` data for a set of segments.

        Parametes
        ---------
        n_iter : int
            Iteration number.
        segments : list of Segment
            Set of segments belonging to iteration `n_iter`.

        """
        segments = sorted(segments, key=attrgetter('seg_id'))
        arrays = [state_to_numpy(s.final_state) for s in segments]
        dtype = arrays[0].dtype  # infer dtype from first segment
        entries = np.fromiter(arrays, dtype=dtype)

        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            n_total_segments = iter_group['seg_index'].shape[0]

            ds = iter_group.require_dataset('final_states', (n_total_segments,), dtype=dtype)
            dsid = ds.id

            msel = h5s.create_simple(entries.shape, (h5s.UNLIMITED,))
            msel.select_all()
            fsel = dsid.get_space()
            for i, segment in enumerate(segments):
                op = h5s.SELECT_OR if i != 0 else h5s.SELECT_SET
                fsel.select_hyperslab((segment.seg_id,), (1,), op=op)

            dsid.write(msel, fsel, entries)

    def write_pcoords(self, n_iter, segments):
        """Write :attr:`~westpa.Segment.pcoord` data for a set of segments.

        Parametes
        ---------
        n_iter : int
            Iteration number.
        segments : list of Segment
            Set of segments belonging to iteration `n_iter`.

        """
        segments = sorted(segments, key=attrgetter('seg_id'))
        # infer shape and dtype from first segment
        pcoord_len, pcoord_ndim = segments[0].pcoord.shape
        pcoord_dtype = segments[0].pcoord.dtype

        with self.lock:
            iter_group = self.get_iter_group(n_iter)
            n_total_segments = iter_group['seg_index'].shape[0]
            n_segments = len(segments)

            default_opts = {'name': 'pcoord', 'h5path': 'pcoord', 'compression': False}
            opts = self.dataset_options.get('pcoord', default_opts)
            shape = (n_total_segments, pcoord_len, pcoord_ndim)
            ds = require_dataset_from_dsopts(iter_group, opts, shape, dtype=pcoord_dtype)
            dsid = ds.id
            entries = np.empty((n_segments, pcoord_len, pcoord_ndim), dtype=pcoord_dtype)

            msel = h5s.create_simple(entries.shape, (h5s.UNLIMITED,) * entries.ndim)
            msel.select_all()
            fsel = dsid.get_space()
            for i, segment in enumerate(segments):
                op = h5s.SELECT_OR if i != 0 else h5s.SELECT_SET
                fsel.select_hyperslab((segment.seg_id, 0, 0), (1, pcoord_len, pcoord_ndim), op=op)
                entries[i] = segment.pcoord

            dsid.write(msel, fsel, entries)

    def get_segments(self, n_iter=None, seg_ids=None, load_pcoords=True, load_auxdata=False):
        segments = super().get_segments(n_iter, seg_ids, load_pcoords=False)
        seg_ids = [s.seg_id for s in segments]

        with self.lock:
            iter_group = self.get_iter_group(n_iter or self.current_iteration)
            if 'initial_states' in iter_group:
                for segment, row in zip(segments, iter_group['initial_states'][seg_ids]):
                    segment.initial_state = state_from_numpy(row)
            if 'final_states' in iter_group:
                for segment, row in zip(segments, iter_group['final_states'][seg_ids]):
                    segment.final_state = state_from_numpy(row)
            if load_pcoords and 'pcoord' in iter_group:
                for segment, row in zip(segments, iter_group['pcoord'][seg_ids]):
                    segment.pcoord = row
            if load_auxdata and 'auxdata' in iter_group:
                for dsname, ds in iter_group['auxdata'].items():
                    for segment, row in zip(segments, ds[seg_ids]):
                        segment.data[dsname] = row

        return segments
