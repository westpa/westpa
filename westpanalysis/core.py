import functools
import operator
import westpa.tools.binning

from abc import ABC, abstractmethod
from westpa.core.h5io import WESTPAH5File
from westpa.core.states import BasisState, InitialState, TargetState


class Run:
    """A completed WESTPA simulation run.

    Parameters
    ----------
    h5filename : str, default 'west.h5'
        Pathname of a WESTPA HDF5 data file.
    name : str, optional
        Name of the run. Default is ''.
    segment_traj_loader : callable, optional
        A function to be used to load segment trajectories. The required
        signature and return type of `seg_traj_loader` are specified by the
        :class:`SegmentTrajectoryLoader` abstract base class.

    """

    DESCRIPTION = 'WESTPA Run'

    def __init__(self, h5filename='west.h5', name=None, segment_traj_loader=None):
        self.h5file = WESTPAH5File(h5filename, 'r')
        self.name = name
        self.segment_traj_loader = segment_traj_loader

    @property
    def segment_traj_loader(self):
        """callable: Function used to load segment trajectories."""
        return self._segment_traj_loader

    @segment_traj_loader.setter
    def segment_traj_loader(self, value):
        if value is not None:
            _ = value(1, 0)
        self._segment_traj_loader = value

    @property
    def _default_name(self):
        return ''

    @property
    def name(self):
        """str: Name of the run."""
        return self._name

    @name.setter
    def name(self, value):
        if value is None:
            self._name = self._default_name
            return
        if not isinstance(value, str):
            raise TypeError('name must be a string')
        self._name = value

    @property
    def summary(self):
        """h5py.Dataset: The 'summary' dataset of the HDF5 file."""
        return self.h5file['summary']

    @property
    def basis_state_info(self):
        """h5py.Dataset: 'bstate_index' dataset."""
        return self.h5file['ibstates']['bstate_index']

    @property
    def num_iterations(self):
        """int: Number of iterations in the run."""
        return len(self.h5file['iterations'])

    @property
    def iterations(self):
        """Iterable[Iteration]: Sequence of iterations."""
        return (Iteration(number, self)
                for number in range(1, self.num_iterations + 1))

    @property
    def num_segments(self):
        """int: Total number of segments."""
        return sum(iteration.num_segments for iteration in self)

    def iteration(self, number):
        """Return the iteration with the given iteration number.

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

    def __container__(self, iteration):
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
        """h5py.Dataset: 'seg_index' dataset of the iteration."""
        return self.h5group['seg_index']

    @property
    def pcoords(self):
        """h5py.Dataset: 'pcoord' dataset of the iteration."""
        return self.h5group['pcoord']

    @property
    def bin_target_counts(self):
        """h5py.Dataset: 'bin_target_counts' dataset of the iteration."""
        return self.h5group.get('bin_target_counts')  # May be None.

    @functools.cached_property
    def bin_mapper(self):
        """BinMapper: Bin mapper used in the iteration."""
        if self.bin_target_counts is None:
            return None
        mapper, _, _ = westpa.tools.binning.mapper_from_hdf5(
            self.run.h5file['bin_topologies'],
            self.h5group.attrs['binhash'])
        return mapper

    @property
    def num_bins(self):
        """int: Number of bins."""
        return self.bin_target_counts.shape[0]

    @property
    def bins(self):
        """Iterable[Bin]: Bins."""
        return (Bin(index, self) for index in range(self.num_bins))

    @property
    def num_segments(self):
        """int: Number of segments in the iteration."""
        return self.segment_info.shape[0]

    @property
    def segment_weights(self):
        """1D ndarray: Statistical weights of the segments."""
        return self.segment_info['weight']

    @property
    def segments(self):
        """Iterable[Segment]: Segments of the iteration."""
        return (Segment(index, self) for index in range(self.num_segments))

    @property
    def _ibstates(self):
        return self.h5group['ibstates']

    @property
    def _tstates(self):
        return self.h5group['tstates']

    @property
    def basis_state_info(self):
        """h5py.Dataset: 'bstate_index' dataset."""
        return self._ibstates['bstate_index']

    @property
    def basis_state_pcoords(self):
        """h5py.Dataset. 'bstate_pcoord' dataset."""
        return self._ibstates['bstate_pcoord']

    @property
    def basis_states(self):
        """list[BasisState]: Basis states in use for the iteration."""
        return [BasisState(info['label'], info['probability'], pcoord=pcoord,
                           auxref=info['auxref'], state_id=state_id)
                for state_id, (info, pcoord) in enumerate(
                    zip(self.basis_state_info, self.basis_state_pcoords))]

    @property
    def initial_state_info(self):
        """h5py.Dataset. 'istate_index' dataset."""
        return self._ibstates['istate_index']

    @property
    def initial_state_pcoords(self):
        """h5py.Dataset. 'istate_pcoord' dataset."""
        return self._ibstates['istate_pcoord']

    @property
    def initial_states(self):
        """list[InitialState]: Initial states."""
        return [InitialState(state_id, info['basis_state_id'],
                             info['iter_created'], iter_used=info['iter_used'],
                             istate_type=info['istate_type'],
                             istate_status=info['istate_status'],
                             pcoord=pcoord)
                for state_id, (info, pcoord) in enumerate(
                    zip(self.initial_state_info, self.initial_state_pcoords))]

    @property
    def target_state_info(self):
        """h5py.Dataset: 'index' dataset for target states."""
        return self._tstates['index']

    @property
    def target_state_pcoords(self):
        """h5py.Dataset: 'pcoord' dataset for target states."""
        return self._tstates['pcoord']

    @property
    def target_states(self):
        """list[TargetState]: Target states."""
        return [TargetState(info['label'], pcoord, state_id=state_id)
                for state_id, (info, pcoord) in enumerate(
                    zip(self.target_state_info, self.target_state_pcoords))]

    @property
    def target(self):
        """Target: Union of bins that serve as the target."""
        if not self.next:
            return None
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

    def segment(self, index):
        """Return the segment with the given index.

        Parameters
        ----------
        index : int
            Segment index (0-based).

        Returns
        -------
        Segment
            The segment indexed by `index`.

        """
        valid_range = range(self.num_segments)
        if index not in valid_range:
            raise ValueError(f'segment index must be in {valid_range}')
        return Segment(index, self)

    def __iter__(self):
        return iter(self.segments)

    def __contains__(self, segment):
        return segment.iteration is self

    def __eq__(self, other):
        return self.number == other.number and self.run == other.run

    def __repr__(self):
        return f'{self.__class__.__name__}({self.number}, {self.run})'


class Segment:
    """A segment of an iteration of a WESTPA simulation.

    Parameters
    ----------
    index : int
        Segment index (0-based).
    iteration : Iteration
        Iteration to which the segment belongs.

    """

    def __init__(self, index, iteration):
        self.index = index
        self.iteration = iteration

    @property
    def run(self):
        """Run: Run to which the segment belongs."""
        return self.iteration.run

    @property
    def info(self):
        """numpy.void: 'seg_index' row of the segment."""
        return self.iteration.segment_info[self.index]

    @property
    def weight(self):
        """float64: Statistical weight of the segment."""
        return self.info['weight']

    @property
    def pcoords(self):
        """2D ndarray: Progress coordinates at each snapshot time."""
        return self.iteration.pcoords[self.index]

    @property
    def num_snapshots(self):
        """int: Number of snapshots (saved frames)."""
        return self.pcoords.shape[0]

    @property
    def parent(self):
        """Segment: The parent of the segment."""
        if not self.iteration.prev:
            return None
        return Segment(self.info['parent_id'], self.iteration.prev)

    @property
    def children(self):
        """Iterable[Segment]: The children of the segment."""
        if not self.iteration.next:
            return ()
        return (segment for segment in self.iteration.next.segments
                if segment.parent == self)

    def trace(self):
        """Return the trace (ancestral line) of the segment.

        Returns
        -------
        Trace
            The trace of the segment.

        """
        return Trace(self)

    def trajectory(self):
        """Return the trajectory of the segment.

        Returns
        -------
        sequence
            The trajectory of the segment.

        See Also
        --------
        Run.segment_traj_loader
            Function used to load the trajectory.

         """
        if not self.run.segment_traj_loader:
            raise ValueError('segment trajectory loader must be set')
        return self.run.segment_traj_loader(
            self.iteration.number, self.index)

    def __len__(self):
        return self.num_snapshots

    def __eq__(self, other):
        return self.index == other.index and self.iteration == other.iteration

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.iteration})'


class SegmentTrajectoryLoader(ABC):
    """API for segment trajectory loaders."""

    @abstractmethod
    def __call__(self, iteration_number, segment_index):
        """Return the trajectory of a given segment.

        Parameters
        ----------
        iteration_number : int
            Number (1-based) of iteration to which the segment belongs.
        segment_index : int
            Index (0-based) of the segment within the iteration.

        Returns
        -------
        sequence
            The trajectory of the segment.

        """
        ...


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

    def __contains__(self, pcoord):
        result = self.mapper.assign([pcoord])
        if result.size != 1:
            raise ValueError('left operand must be a single point in '
                             'progress coordinate space')
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
    """A trace of a segment back to a source or initial state.

    Parameters
    ----------
    segment : Segment
        The traced segment.
    source : Bin or BinUnion, optional
        The source (macro)state. If provided, the trace is continued
        only as far back as `source`. Otherwise, the trace extends back
        to the initial state.

    """

    def __init__(self, segment, source=None):
        if source is None:
            source = BinUnion()

        segments = []
        while segment and segment.pcoords[-1] not in source:
            segments.append(segment)
            segment = segment.parent

        self.segments = tuple(reversed(segments))
        self.source = source

    def trajectory(self, concat_operator=None):
        """Return the trajectory of the trace.

        Parameters
        ----------
        concat_operator : callable, default :func:`operator.concat`
            A function that takes two trajectories and returns their
            concatenation. This parameter may be provided if the type
            returned by the segment trajectory loader does not support
            concatenation via ``traj1 + traj2`` (the default).

        Returns
        -------
        sequence
            The trajectory of the trace.

        See Also
        --------
        Run.segment_traj_loader
            Function used to load the segment trajectories.

        """
        if concat_operator is None:
            concat_operator = operator.concat
        trajectories = (segment.trajectory() for segment in self.segments)
        return functools.reduce(concat_operator, trajectories)

    def __repr__(self):
        s = f'{self.__class__.__name__}({self.segments[-1]}'
        if self.source:
            s += f', source={self.source}'
        return s + ')'
