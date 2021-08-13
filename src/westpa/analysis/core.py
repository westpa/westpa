import itertools
import pandas as pd
import sys
import westpa.tools.binning

from functools import cached_property
from westpa.core.h5io import WESTPAH5File
from westpa.core.states import BasisState, InitialState, TargetState


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
    def successful_walkers(self):
        """Iterable[Walker]: Walkers that stopped in the target."""
        return itertools.chain.from_iterable(iteration.successful_walkers for iteration in self)

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
        if self.number == self.run.num_iterations:
            return None
        return self.run.iteration(self.number + 1)

    @property
    def segment_info(self):
        """pd.DataFrame: Segment summary data for the iteration."""
        return pd.DataFrame(self.h5group['seg_index'][:], dtype=object)

    @cached_property
    def pcoords(self):
        """3D ndarray: Progress coordinate snaphots of each walker."""
        return self.h5group['pcoord'][:]

    @cached_property
    def weights(self):
        """1D ndarray: Statistical weight of each walker."""
        return self.segment_info['weight']

    @property
    def bin_target_counts(self):
        """1D ndarray, dtype=uint64: Target count for each bin."""
        val = self.h5group.get('bin_target_counts')
        if val is None:
            return None
        return val[:]

    @cached_property
    def bin_mapper(self):
        """BinMapper: Bin mapper used in the iteration."""
        if self.bin_target_counts is None:
            return None
        mapper, _, _ = westpa.tools.binning.mapper_from_hdf5(self.run.h5file['bin_topologies'], self.h5group.attrs['binhash'])
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
        return (Bin(index, self) for index in range(self.num_bins))

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
    def successful_walkers(self):
        """Iterable[Walker]: Walkers that stopped in the target."""
        return (walker for walker in self if walker.successful)

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

    @cached_property
    def target(self):
        """Target: Union of bins that serve as the target."""
        return Target(self)

    def bin(self, index):
        """Return the bin with the given index.

        Parameters
        ----------
        index : int
            Bin index.

        Returns
        -------
        Bin
            The bin indexed by `index`.

        """
        valid_range = range(self.num_bins)
        if index not in valid_range:
            raise ValueError(f'bin index must be in {valid_range}')
        return Bin(index, self)

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
        return self.iteration.pcoords[self.index]

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
        return Walker(self.segment_info['parent_id'], self.iteration.prev)

    @property
    def children(self):
        """Iterable[Walker]: The children of the walker."""
        if not self.iteration.next:
            return ()
        return (walker for walker in self.iteration.next.walkers if walker.parent == self)

    @property
    def successful(self):
        """bool: True if the walker stopped in the target, False otherwise."""
        return self.stopped_in(self.iteration.target)

    def stopped_in(self, pcoord_subset):
        """Return True if the walker stopped (i.e., terminated) in a given
        subset of progress coordinate space, False otherwise.

        Parameters
        ----------
        pcoord_subset : Bin, BinUnion, or Container
            A :type:`Container` object representing a subset of progress
            coordinate space.

        Returns
        -------
        bool
            Whether the walker stopped in `pcoord_subset`.

        """
        return self.pcoords[-1] in pcoord_subset

    def trace(self, source=None):
        """Return the trace (ancestral line) of the walker.

        For full documentation see :class:`Trace`.

        Returns
        -------
        Trace
            The trace of the walker.

        """
        return Trace(self, source=source)

    def __eq__(self, other):
        return self.index == other.index and self.iteration == other.iteration

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.iteration})'


class Bin:
    """A bin used in the resampling step of an iteration.

    Parameters
    ----------
    index : int
        Bin index.
    iteration : Iteration
        The iteration in which the bin was used.

    """

    def __init__(self, index, iteration):
        self.index = index
        self.iteration = iteration

    @property
    def mapper(self):
        """BinMapper: Bin mapper that defines the bin."""
        return self.iteration.bin_mapper

    @property
    def target_count(self):
        """int: Target number of particles in the bin."""
        return self.iteration.bin_target_counts[self.index]

    def union(self, *others):
        """Return the union of this bin with other bins.

        Parameters
        ----------
        *others : Bin
            Other bins comprising the union.

        Return
        ------
        BinUnion
            The union of `self` and `others`.

        """
        return BinUnion(self, *others)

    def __or__(self, other):
        return self.union(other)

    def __contains__(self, pcoord):
        result = self.mapper.assign([pcoord])
        if result.size != 1:
            raise ValueError('left operand must be a single point in progress coordinate space')
        return result[0] == self.index

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.iteration})'


class BinUnion:
    """A union of bins.

    Parameters
    ----------
    *bins : Bin
        The bins comprising the union.

    """

    def __init__(self, *bins):
        if not all(isinstance(bin_, Bin) for bin_ in bins):
            raise TypeError('arguments must be of type Bin')
        self._bins = bins

    @property
    def bins(self):
        """tuple[Bin]: Bins comprising the union."""
        return self._bins

    def union(self, *others):
        """Return the union of this bin with other bins.

        Parameters
        ----------
        *others : BinUnion or Bin
            Other bins comprising the union.

        Return
        ------
        BinUnion
            The union of `self` and `others`.

        """
        bins = list(self.bins)
        for other in others:
            if isinstance(other, Bin):
                bins.append(other)
                continue
            bins += other.bins
        return BinUnion(*bins)

    def __or__(self, other):
        return self.union(other)

    def __bool__(self):
        return bool(self.bins)

    def __contains__(self, pcoord):
        return any(pcoord in bin_ for bin_ in self.bins)

    def __repr__(self):
        return f'{self.__class__.__name__}{self.bins}'


class Target(BinUnion):
    """A union of bins that serve as the target for an iteration.

    Parameters
    ----------
    iteration : Iteration
        The iteration for which the target is defined.

    Notes
    -----
    The bins of the target belong to the *next* iteration.

    """

    def __init__(self, iteration):
        if not iteration.next or not iteration.target_states:
            super().__init__()
        else:
            pcoords = iteration.target_state_pcoords[:]
            bin_indices = set(iteration.next.bin_mapper.assign(pcoords))
            super().__init__(*(Bin(i, iteration.next) for i in bin_indices))
        self.iteration = iteration

    @property
    def states(self):
        """list[TargetState]: Target states defining the target."""
        return self.iteration.target_states

    def __repr__(self):
        return f'{self.__class__.__name__}({self.iteration})'


class Trace:
    """A trace of a walker's ancestry back to a source or initial state.

    Parameters
    ----------
    walker : Walker
        The terminal walker.
    source : Bin, BinUnion, or Container, optional
        The source (macro)state, specified as a :type:`Container` object
        whose :meth:`__contains__` method is the indicator function for
        the corresponding subset of progress coordinate space. If `source`
        is provided, the trace is continued only as far back as the last
        walker that stopped in `source`. Otherwise, the trace extends back
        to the initial state.

    """

    def __init__(self, walker, source=None, max_length=None):
        if source is None:
            source = BinUnion()

        if max_length is None:
            max_length = sys.maxsize
        else:
            if max_length <= 0:
                raise ValueError('max_length must be greater than 0')
            if not isinstance(max_length, int):
                raise TypeError('max_length must be an integer')

        walkers = []
        while walker and walker.pcoords[-1] not in source:
            walkers.append(walker)
            if len(walkers) == max_length:
                break
            walker = walker.parent

        self.walkers = tuple(reversed(walkers))
        self.source = source
        self.max_length = max_length

    def __len__(self):
        return len(self.walkers)

    def __iter__(self):
        return iter(self.walkers)

    def __contains__(self, walker):
        return walker in self.walkers

    def __repr__(self):
        s = f'Trace({self.walkers[-1]}'
        if self.source:
            s += f',\n      source={self.source}'
        if self.max_length < sys.maxsize:
            s += f',\n      max_length={self.max_length}'
        return s + ')'
