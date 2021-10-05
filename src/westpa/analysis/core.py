import itertools
import sys

import numpy as np
import pandas as pd

from westpa.core.binning.assign import BinMapper
from westpa.core.h5io import WESTPAH5File
from westpa.core.segment import Segment
from westpa.core.states import BasisState, InitialState, TargetState

from westpa.tools.binning import mapper_from_hdf5


class Run:
    """A completed WESTPA simulation run.

    Parameters
    ----------
    h5filename : str or file-like object, default 'west.h5'
        Pathname or stream of a WESTPA HDF5 data file.
    name : str, optional
        Name of the run. Default is the empty string.

    """

    DESCRIPTION = 'WESTPA Run'

    def __init__(self, h5filename='west.h5', name=None):
        self.h5filename = h5filename
        self.name = name

    @property
    def h5filename(self):
        return self.h5file.filename

    @h5filename.setter
    def h5filename(self, value):
        self.h5file = WESTPAH5File(value, 'r')

    def __getstate__(self):
        state = self.__dict__.copy()
        state['h5filename'] = state['h5file'].filename
        del state['h5file']
        return state

    def __setstate__(self, state):
        state['h5file'] = WESTPAH5File(state['h5filename'], 'r')
        del state['h5filename']
        self.__dict__.update(state)

    @property
    def name(self):
        """str: Name of the run."""
        return self._name

    @name.setter
    def name(self, value):
        if value is None:
            value = ''
        elif not isinstance(value, str):
            raise TypeError('name must be a string')
        self._name = value

    @property
    def summary(self):
        """pd.DataFrame: Summary data for the run."""
        return pd.DataFrame(self.h5file['summary'][:])

    @property
    def num_iterations(self):
        """int: Number of completed iterations."""
        return self.h5file.attrs['west_current_iteration'] - 1

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
        """int: Alias self.num_walkers."""
        return self.num_walkers

    @property
    def walkers(self):
        """Iterable[Walker]: All walkers in the run."""
        return itertools.chain.from_iterable(self)

    @property
    def recycled_walkers(self):
        """Iterable[Walker]: Walkers that stopped in the sink."""
        return itertools.chain.from_iterable(iteration.recycled_walkers for iteration in self)

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
        return iteration.run is self

    def __eq__(self, other):
        return self.h5file == other.h5file

    def __repr__(self):
        if self.name:
            return f'<{self.DESCRIPTION} "{self.name}" at {hex(id(self))}>'
        return f'<{self.DESCRIPTION} at {hex(id(self))}>'


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
        return self.run.iteration(self.number + 1)

    @property
    def segment_info(self):
        """pd.DataFrame: Segment summary data for the iteration."""
        return pd.DataFrame(self.h5group['seg_index'][:], dtype=object)

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
        return self.segment_info.shape[0]

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
    def _ibstates(self):
        return self.h5group['ibstates']

    @property
    def _tstates(self):
        return self.h5group.get('tstates')  # None for equilibrium sampling

    @property
    def basis_state_info(self):
        """pd.DataFrame: Basis state summary data."""
        return pd.DataFrame(self._ibstates['bstate_index'][:])

    @property
    def basis_state_pcoords(self):
        """2D ndarray: Progress coordinates of each basis state."""
        return self._ibstates['bstate_pcoord'][:]

    @property
    def basis_states(self):
        """list[BasisState]: Basis states in use for the iteration."""
        return [
            BasisState(row.label, row.probability, pcoord=pcoord, auxref=row.auxref, state_id=row.Index)
            for row, pcoord in zip(self.basis_state_info.itertuples(), self.basis_state_pcoords)
        ]

    @property
    def initial_state_info(self):
        """pd.DataFrame: Initial state summary data."""
        return pd.DataFrame(self._ibstates['istate_index'][:])

    @property
    def initial_state_pcoords(self):
        """2D ndarray: Progress coordinates of each initial state."""
        return self._ibstates['istate_pcoord']

    @property
    def initial_states(self):
        """list[InitialState]: Initial states."""
        return [
            InitialState(
                row.Index,
                row.basis_state_id,
                row.iter_created,
                iter_used=row.iter_used,
                istate_type=row.istate_type,
                istate_status=row.istate_status,
                pcoord=pcoord,
            )
            for row, pcoord in zip(self.initial_state_info.itertuples(), self.initial_state_pcoords)
        ]

    @property
    def target_state_info(self):
        """pd.DataFrame: Target state summary data."""
        return pd.DataFrame(self._tstates['index'][:]) if self._tstates else None

    @property
    def target_state_pcoords(self):
        """2D ndarray: Progress coordinates of each target state."""
        return self._tstates['pcoord'][:] if self._tstates else None

    @property
    def target_states(self):
        """list[TargetState]: Target states."""
        if self._tstates is None:
            return []
        return [
            TargetState(row.label, pcoord, state_id=row.Index)
            for row, pcoord in zip(self.target_state_info.itertuples(), self.target_state_pcoords)
        ]

    @property
    def sink(self):
        """BinUnion or None: Union of bins serving as the recycling sink."""
        if not self.target_states:
            return None
        mapper = self.next.bin_mapper
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

    def __iter__(self):
        return iter(self.walkers)

    def __contains__(self, walker):
        return walker.iteration is self

    def __eq__(self, other):
        return self.number == other.number and self.run == other.run

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
    def segment_info(self):
        """pd.Series: Segment summary data."""
        return self.iteration.segment_info.iloc[self.index]

    @property
    def parent(self):
        """Walker: The parent of the walker."""
        if not self.iteration.prev:
            return None
        index = self.iteration.h5group['seg_index']['parent_id'][self.index]
        return Walker(index, self.iteration.prev)

    @property
    def children(self):
        """Iterable[Walker]: The children of the walker."""
        if not self.iteration.next:
            return ()
        return (walker for walker in self.iteration.next.walkers if walker.parent == self)

    @property
    def recycled(self):
        """bool: True if the walker stopped in the sink, False otherwise."""
        endpoint_type = self.iteration.h5group['seg_index']['endpoint_type'][self.index]
        return endpoint_type == Segment.SEG_ENDPOINT_RECYCLED

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
            Other :type:`BinUnion` instances, consisting of bins defined by
            the same underlying bin mapper.

        Return
        ------
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
            Other :type:`BinUnion` instances, consisting of bins defined by
            the same underlying bin mapper.

        Return
        ------
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
    """A trace of a walker's ancestry back to a source or initial state.

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
        if source is None:
            source = ()

        if max_length is None:
            max_length = sys.maxsize
        else:
            max_length = int(max_length)
            if max_length < 1:
                raise ValueError('max_length must be at least 1')

        walkers = []
        while walker and (walker.pcoords[-1] not in source) and (len(walkers) < max_length):
            walkers.append(walker)
            walker = walker.parent

        self.walkers = tuple(reversed(walkers))
        self.source = source
        self.max_length = max_length

    @property
    def pcoords(self):
        """2D ndarray: Progress coordinate snapshots."""
        return np.concatenate([walker.pcoords for walker in self])

    @property
    def num_snapshots(self):
        """int: Number of snapshots."""
        return sum(walker.num_snapshots for walker in self)

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
            s += f',\n      source={self.source}'
        if self.max_length < sys.maxsize:
            s += f',\n      max_length={self.max_length}'
        return s + ')'
