import itertools
import numpy as np
import pandas as pd
import sys

from westpa.core.binning.assign import BinMapper
from westpa.core.h5io import WESTPAH5File, tostr
from westpa.core.segment import Segment
from westpa.core.states import BasisState, InitialState, TargetState
from westpa.tools.binning import mapper_from_hdf5


class Run:
    """A read-only view of a WESTPA simulation run.

    Parameters
    ----------
    h5filename : str or file-like object, default 'west.h5'
        Pathname or stream of a main WESTPA HDF5 data file.

    """

    DESCRIPTION = 'WESTPA Run'

    def __init__(self, h5filename='west.h5'):
        self.h5filename = h5filename

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        self.close()

    @classmethod
    def open(cls, h5filename='west.h5'):
        """Alternate constructor.

        Parameters
        ----------
        h5filename : str or file-like object, default 'west.h5'
            Pathname or stream of a main WESTPA HDF5 data file.

        """
        return cls(h5filename)

    def close(self):
        """Close the Run instance by closing the underlying WESTPA HDF5 file."""
        self.h5file.close()

    @property
    def closed(self):
        """bool: Whether the Run instance is closed."""
        return not bool(self.h5file)

    @property
    def h5filename(self):
        return self.h5file.filename

    @h5filename.setter
    def h5filename(self, value):
        try:
            h5file = WESTPAH5File(value, 'r')
        except FileNotFoundError as e:
            e.strerror = f'Failed to open {self.DESCRIPTION}: file {value!r} not found'
            raise e.with_traceback(None)
        self.h5file = h5file

    @property
    def summary(self):
        """pd.DataFrame: Summary data by iteration."""
        df = pd.DataFrame(
            self.h5file['summary'][: self.num_iterations],
            index=range(1, self.num_iterations + 1),
            dtype=object,
        )
        df.pop('norm')  # should always be 1.0
        df.pop('binhash')  # not human readable
        return df

    @property
    def num_iterations(self):
        """int: Number of completed iterations."""
        if not hasattr(self, '_num_iterations'):
            current = self.h5file.attrs['west_current_iteration']
            grp = self.h5file.get_iter_group(current)
            if (grp['seg_index']['status'] == Segment.SEG_STATUS_COMPLETE).all():
                self._num_iterations = current
            else:
                self._num_iterations = current - 1
        return self._num_iterations

    @property
    def iterations(self):
        """Sequence[Iteration]: Sequence of iterations."""
        return [Iteration(number, self) for number in range(1, self.num_iterations + 1)]

    @property
    def num_walkers(self):
        """int: Total number of walkers."""
        return sum(iteration.num_segments for iteration in self)

    @property
    def num_segments(self):
        """int: Total number of trajectory segments (alias self.num_walkers)."""
        return self.num_walkers

    @property
    def walkers(self):
        """Iterable[Walker]: All walkers in the run."""
        return itertools.chain.from_iterable(self)

    @property
    def recycled_walkers(self):
        """Iterable[Walker]: Walkers that stopped in the sink."""
        return itertools.chain.from_iterable(iteration.recycled_walkers for iteration in self)

    @property
    def initial_walkers(self):
        """Iterable[Walker]: Walkers whose parents are initial states."""
        return itertools.chain.from_iterable(iteration.initial_walkers for iteration in self)

    def iteration(self, number):
        """Return a specific iteration.

        Parameters
        ----------
        number : int
            Iteration number (1-based).

        Returns
        -------
        Iteration
            The iteration indexed by `number`.

        """
        valid_range = range(1, self.num_iterations + 1)
        if number not in valid_range:
            raise ValueError(f'iteration number must be in {valid_range}')
        return Iteration(number, self)

    def __len__(self):
        return self.num_iterations

    def __iter__(self):
        return iter(self.iterations)

    def __contains__(self, iteration):
        return iteration.run == self

    def __eq__(self, other):
        return self.h5file == other.h5file

    def __hash__(self):
        return hash(self.h5file)

    def __bool__(self):
        return not self.closed

    def __repr__(self):
        if self.closed:
            return f'<Closed {self.DESCRIPTION} at {hex(id(self))}>'
        return f'<{self.DESCRIPTION} with {self.num_iterations} iterations at {hex(id(self))}>'


class Iteration:
    """An iteration of a WESTPA simulation.

    Parameters
    ----------
    number : int
        Iteration number (1-based).
    run : Run
        Simulation run to which the iteration belongs.

    """

    def __init__(self, number, run):
        self.number = number
        self.run = run

    @property
    def h5group(self):
        """h5py.Group: HDF5 group containing the iteration data."""
        return self.run.h5file.get_iter_group(self.number)

    @property
    def prev(self):
        """Iteration: Previous iteration."""
        if self.number == 1:
            return None
        return self.run.iteration(self.number - 1)

    @property
    def next(self):
        """Iteration: Next iteration."""
        if self.number == self.run.num_iterations:
            return None
        return self.run.iteration(self.number + 1)

    @property
    def summary(self):
        """pd.DataFrame: Iteration summary."""
        df = pd.DataFrame(
            self.run.h5file['summary'][[self.number - 1]],
            index=[self.number],
            dtype=object,
        )
        df.pop('norm')  # should always be 1.0
        df.pop('binhash')  # not human readable
        return df.iloc[0]

    @property
    def segment_summaries(self):
        """pd.DataFrame: Segment summary data for the iteration."""
        df = pd.DataFrame(self.h5group['seg_index'][:], dtype=object)

        # Make 'endpoint_type' and 'status' human-readable.
        names = map(Segment.endpoint_type_names.get, df['endpoint_type'])
        df['endpoint_type'] = [name.split('_')[-1] for name in names]
        names = map(Segment.status_names.get, df['status'])
        df['status'] = [name.split('_')[-1] for name in names]

        return df

    @property
    def pcoords(self):
        """3D ndarray: Progress coordinate snaphots of each walker."""
        return self.h5group['pcoord'][:]

    @property
    def weights(self):
        """1D ndarray: Statistical weight of each walker."""
        return self.h5group['seg_index']['weight']

    @property
    def bin_target_counts(self):
        """1D ndarray, dtype=uint64: Target count for each bin."""
        val = self.h5group.get('bin_target_counts')
        if val is None:
            return None
        return val[:]

    @property
    def bin_mapper(self):
        """BinMapper: Bin mapper used in the iteration."""
        if self.bin_target_counts is None:
            return None
        mapper, _, _ = mapper_from_hdf5(self.run.h5file['bin_topologies'], self.h5group.attrs['binhash'])
        return mapper

    @property
    def num_bins(self):
        """int: Number of bins."""
        if self.number == 1:
            return 1
        return self.bin_target_counts.shape[0]

    @property
    def bins(self):
        """Iterable[Bin]: Bins."""
        mapper = self.bin_mapper
        return (Bin(index, mapper) for index in range(self.num_bins))

    @property
    def num_walkers(self):
        """int: Number of walkers in the iteration."""
        return self.h5group['seg_index'].shape[0]

    @property
    def num_segments(self):
        """int: Number of trajectory segments (alias self.num_walkers)."""
        return self.num_walkers

    @property
    def walkers(self):
        """Iterable[Walker]: Walkers in the iteration."""
        return (Walker(index, self) for index in range(self.num_walkers))

    @property
    def recycled_walkers(self):
        """Iterable[Walker]: Walkers that stopped in the sink."""
        endpoint_type = self.h5group['seg_index']['endpoint_type']
        indices = np.flatnonzero(endpoint_type == Segment.SEG_ENDPOINT_RECYCLED)
        return (Walker(index, self) for index in indices)

    @property
    def initial_walkers(self):
        """Iterable[Walker]: Walkers whose parents are initial states."""
        parent_ids = self.h5group['seg_index']['parent_id']
        return (walker for walker, parent_id in zip(self, parent_ids) if parent_id < 0)

    @property
    def auxiliary_data(self):
        """h5py.Group or None: Auxiliary data stored for the iteration."""
        return self.h5group.get('auxdata')

    @property
    def basis_state_summaries(self):
        """pd.DataFrame: Basis state summary data."""
        df = pd.DataFrame(self.h5group['ibstates']['bstate_index'][:])
        df['label'] = tostr(df['label'].str)
        df['auxref'] = tostr(df['auxref'].str)
        return df

    @property
    def basis_state_pcoords(self):
        """2D ndarray: Progress coordinates of each basis state."""
        return self.h5group['ibstates']['bstate_pcoord'][:]

    @property
    def basis_states(self):
        """list[BasisState]: Basis states in use for the iteration."""
        return [
            BasisState(row.label, row.probability, pcoord=pcoord, auxref=row.auxref, state_id=row.Index)
            for row, pcoord in zip(self.basis_state_summaries.itertuples(), self.basis_state_pcoords)
        ]

    @property
    def has_target_states(self):
        """bool: Whether target (sink) states are defined for this iteration."""
        return 'tstates' in self.h5group

    @property
    def target_state_summaries(self):
        """pd.DataFrame or None: Target state summary data."""
        if self.has_target_states:
            df = pd.DataFrame(self.h5group['tstates']['index'][:])
            df['label'] = tostr(df['label'].str)

            return df
        else:
            return None

    @property
    def target_state_pcoords(self):
        """2D ndarray or None: Progress coordinates of each target state."""
        return self.h5group['tstates']['pcoord'][:] if self.has_target_states else None

    @property
    def target_states(self):
        """list[TargetState]: Target states in use for the iteration."""
        if not self.has_target_states:
            return []
        return [
            TargetState(row.label, pcoord, state_id=row.Index)
            for row, pcoord in zip(self.target_state_summaries.itertuples(), self.target_state_pcoords)
        ]

    @property
    def sink(self):
        """BinUnion or None: Union of bins serving as the recycling sink."""
        if not self.has_target_states:
            return None
        mapper = Iteration(self.number + 1, self.run).bin_mapper
        return BinUnion(mapper.assign(self.target_state_pcoords), mapper)

    def bin(self, index):
        """Return the bin with the given index.

        Parameters
        ----------
        index : int
            Bin index (0-based).

        Returns
        -------
        Bin
            The bin indexed by `index`.

        """
        return Bin(index, self.bin_mapper)

    def walker(self, index):
        """Return the walker with the given index.

        Parameters
        ----------
        index : int
            Walker index (0-based).

        Returns
        -------
        Walker
            The walker indexed by `index`.

        """
        valid_range = range(self.num_walkers)
        if index not in valid_range:
            raise ValueError(f'walker index must be in {valid_range}')
        return Walker(index, self)

    def basis_state(self, index):
        """Return the basis state with the given index.

        Parameters
        ----------
        index : int
            Basis state index (0-based).

        Returns
        -------
        BasisState
            The basis state indexed by `index`.

        """
        row = self.h5group['ibstates']['bstate_index'][index]
        pcoord = self.h5group['ibstates']['bstate_pcoord'][index]
        return BasisState(
            tostr(row['label']),
            row['probability'],
            pcoord=pcoord,
            auxref=tostr(row['auxref']),
            state_id=index,
        )

    def target_state(self, index):
        """Return the target state with the given index.

        Parameters
        ----------
        index : int
            Target state index (0-based).

        Returns
        -------
        TargetState
            The target state indexed by `index`.

        """
        row = self.h5group['tstates']['index'][index]
        pcoord = self.h5group['tstates']['pcoord'][index]
        return TargetState(tostr(row['label']), pcoord, state_id=index)

    def __iter__(self):
        return iter(self.walkers)

    def __contains__(self, walker):
        return walker.iteration == self

    def __eq__(self, other):
        return self.number == other.number and self.run == other.run

    def __hash__(self):
        return hash((self.number, self.run))

    def __repr__(self):
        return f'{self.__class__.__name__}({self.number}, {self.run})'


class Walker:
    """A walker in an iteration of a WESTPA simulation.

    Parameters
    ----------
    index : int
        Walker index (0-based).
    iteration : Iteration
        Iteration to which the walker belongs.

    """

    def __init__(self, index, iteration):
        self.index = index
        self.iteration = iteration

    @property
    def run(self):
        """Run: Run to which the walker belongs."""
        return self.iteration.run

    @property
    def weight(self):
        """float64: Statistical weight of the walker."""
        return self.iteration.weights[self.index]

    @property
    def pcoords(self):
        """2D ndarray: Progress coordinate snapshots."""
        return self.iteration.h5group['pcoord'][self.index]

    @property
    def num_snapshots(self):
        """int: Number of snapshots."""
        return self.pcoords.shape[0]

    @property
    def segment_summary(self):
        """pd.Series: Segment summary data."""
        df = pd.DataFrame(
            self.iteration.h5group['seg_index'][[self.index]],
            index=[self.index],
            dtype=object,
        )

        # Make 'endpoint_type' and 'status' human-readable.
        names = map(Segment.endpoint_type_names.get, df['endpoint_type'])
        df['endpoint_type'] = [name.split('_')[-1] for name in names]
        names = map(Segment.status_names.get, df['status'])
        df['status'] = [name.split('_')[-1] for name in names]

        return df.iloc[0]

    @property
    def parent(self):
        """Walker or InitialState: The parent of the walker."""
        parent_id = self.iteration.h5group['seg_index']['parent_id'][self.index]

        if parent_id >= 0:
            return Walker(parent_id, self.iteration.prev)

        istate_id = -(parent_id + 1)
        row = self.iteration.h5group['ibstates']['istate_index'][istate_id]

        # Initial states may or may not be generated from a basis state.
        bstate_id = row['basis_state_id']
        try:
            bstate = self.iteration.basis_state(bstate_id)
        except IndexError:
            bstate = None
            bstate_id = None

        return InitialState(
            istate_id,
            bstate_id,
            row['iter_created'],
            iter_used=row['iter_used'],
            istate_type=row['istate_type'],
            istate_status=row['istate_status'],
            pcoord=self.iteration.h5group['ibstates']['istate_pcoord'][istate_id],
            basis_state=bstate,
        )

    @property
    def children(self):
        """Iterable[Walker]: The children of the walker."""
        next = self.iteration.next
        if next is None:
            return ()
        indices = np.flatnonzero(next.h5group['seg_index']['parent_id'] == self.index)
        return (Walker(index, next) for index in indices)

    @property
    def recycled(self):
        """bool: True if the walker stopped in the sink, False otherwise."""
        endpoint_type = self.iteration.h5group['seg_index']['endpoint_type'][self.index]
        return endpoint_type == Segment.SEG_ENDPOINT_RECYCLED

    @property
    def initial(self):
        """bool: True if the parent of the walker is an initial state, False otherwise."""
        return self.iteration.h5group['seg_index']['parent_id'][self.index] < 0

    @property
    def auxiliary_data(self):
        """dict: Auxiliary data for the walker."""
        data = self.iteration.auxiliary_data or {}
        return {name: data[name][self.index] for name in data}

    def trace(self, **kwargs):
        """Return the trace (ancestral line) of the walker.

        For full documentation see :class:`Trace`.

        Returns
        -------
        Trace
            The trace of the walker.

        """
        return Trace(self, **kwargs)

    def __eq__(self, other):
        return self.index == other.index and self.iteration == other.iteration

    def __hash__(self):
        return hash((self.index, self.iteration))

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.iteration})'


class BinUnion:
    """A (disjoint) union of bins defined by a common bin mapper.

    Parameters
    ----------
    indices : iterable of int
        The indices of the bins comprising the union.
    mapper : BinMapper
        The bin mapper defining the bins.

    """

    def __init__(self, indices, mapper):
        if not isinstance(mapper, BinMapper):
            raise TypeError(f'mapper must be an instance of {BinMapper}')

        indices = set(indices)
        valid_range = range(mapper.nbins)
        if any(index not in valid_range for index in indices):
            raise ValueError(f'bin indices must be in {valid_range}')

        self.indices = indices
        self.mapper = mapper

    def union(self, *others):
        """Return the union of the bin union and all others.

        Parameters
        ----------
        *others : BinUnion
            Other :class:`BinUnion` instances, consisting of bins defined by
            the same underlying bin mapper.

        Returns
        -------
        BinUnion
            The union of `self` and `others`.

        """
        if any(other.mapper != self.mapper for other in others):
            raise ValueError('bins must be defined by the same bin mapper')
        indices = self.indices.union(*(other.indices for other in others))
        return BinUnion(indices, self.mapper)

    def intersection(self, *others):
        """Return the intersection of the bin union and all others.

        Parameters
        ----------
        *others : BinUnion
            Other :class:`BinUnion` instances, consisting of bins defined by
            the same underlying bin mapper.

        Returns
        -------
        BinUnion
            The itersection of `self` and `others`.

        """
        if any(other.mapper != self.mapper for other in others):
            raise ValueError('bins must be defined by the same bin mapper')
        indices = self.indices.intersection(*(other.indices for other in others))
        return BinUnion(indices, self.mapper)

    def __contains__(self, coord):
        result = self.mapper.assign([coord])
        if result.size != 1:
            raise ValueError('left operand must be a single coordinate tuple')
        return result[0] in self.indices

    def __or__(self, other):
        return self.union(other)

    def __and__(self, other):
        return self.intersection(other)

    def __bool__(self):
        return bool(self.indices)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.indices}, {self.mapper})'


class Bin(BinUnion):
    """A bin defined by a bin mapper.

    Parameters
    ----------
    index : int
        The index of the bin.
    mapper : BinMapper
        The bin mapper defining the bin.

    """

    def __init__(self, index, mapper):
        super().__init__({index}, mapper)
        self.index = index

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.mapper})'


class Trace:
    """A trace of a walker's ancestry.

    Parameters
    ----------
    walker : Walker
        The terminal walker.
    source : Bin, BinUnion, or collections.abc.Container, optional
        A source (macro)state, specified as a container object whose
        :meth:`__contains__` method is the indicator function for the
        corresponding subset of progress coordinate space. The trace is
        stopped upon encountering a walker that stopped in `source`.
    max_length : int, optional
        The maximum number of walkers in the trace.

    """

    def __init__(self, walker, source=None, max_length=None):
        if max_length is None:
            max_length = sys.maxsize
        else:
            max_length = int(max_length)
            if max_length < 1:
                raise ValueError('max_length must be at least 1')

        walkers = []
        initial_state = None
        while len(walkers) < max_length:
            if source and walker.pcoords[-1] in source:
                break
            walkers.append(walker)
            parent = walker.parent
            if isinstance(parent, InitialState):
                initial_state = parent
                break
            walker = parent
        walkers.reverse()

        self.walkers = walkers
        self.initial_state = initial_state
        self.source = source
        self.max_length = max_length

    def __len__(self):
        return len(self.walkers)

    def __iter__(self):
        return iter(self.walkers)

    def __contains__(self, walker):
        return walker in self.walkers

    def __getitem__(self, key):
        return self.walkers[key]

    def __repr__(self):
        s = f'Trace({self.walkers[-1]}'
        if self.source:
            s += f', source={self.source}'
        if self.max_length < sys.maxsize:
            s += f', max_length={self.max_length}'
        return s + ')'
