import functools
import operator

from abc import ABC, abstractmethod
from westpa.core.h5io import WESTPAH5File


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
    def dataset(self):
        """WESTDataset: Simulation record to which the segment belongs."""
        return self.iteration.dataset

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
        """2D ndarray: Progress coordinate values of each snapshot."""
        return self.iteration.pcoord[self.index]

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
            return []
        return [segment for segment in self.iteration.next.segments
                if segment.parent == self]

    def trace(self):
        """Return the trace (ancestral line) of the segment.

        Returns
        -------
        Iterable[Segment]

        See :func:`trace`.

        """
        return trace(self)

    def trajectory(self):
        """Return the trajectory of the segment.

        Returns
        -------
        sequence
            The trajectory of the segment.

        See Also
        --------
        :attr:`WESTDataset.segment_traj_loader`
            Function used to load the trajectory.

         """
        if not self.dataset.segment_traj_loader:
            return None
        return self.dataset.segment_traj_loader(
            self.iteration.number, self.index)

    def trace_trajectory(self, concat_operator=None):
        """Return the trajectory of the trace of the segment.

        Parameters
        ----------
        concat_operator : callable, default :func:`operator.concat`
            A function that takes two trajectories and returns their
            concatenation. This parameter may be provided if the type
            returned by :func:`self.dataset.seg_traj_loader` does not
            support concatenation via `traj1 + traj2` (the default).

        Returns
        -------
        sequence
            The trajectory of the trace (ancestral line) of the segment.

        See Also
        --------
        :attr:`WESTDataset.segment_traj_loader`
            Function used to load the segment trajectories.

        """
        if not self.dataset.segment_traj_loader:
            return None
        if concat_operator is None:
            concat_operator = operator.concat
        trajectories = (segment.trajectory() for segment in self.trace())
        return functools.reduce(concat_operator, trajectories)

    def __len__(self):
        return self.num_snapshots

    def __eq__(self, other):
        return self.index == other.index and self.iteration == other.iteration

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.iteration})'


class Iteration:
    """An iteration of a WESTPA simulation.

    Parameters
    ----------
    number : int
        Iteration number (1-based).
    dataset : WESTDataset
        Simulation record to which the iteration belongs.

    """
    def __init__(self, number, dataset):
        self.number = number
        self.dataset = dataset

    @property
    def group(self):
        """h5py.Group: HDF5 group containing the iteration data."""
        return self.dataset.iteration_group(self.number)

    @property
    def segment_info(self):
        """h5py.Dataset: 'seg_index' dataset of the iteration."""
        return self.group['seg_index']

    @property
    def pcoord(self):
        """h5py.Dataset: 'pcoord' dataset of the iteration."""
        return self.group['pcoord']

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

    def segment(self, index):
        """Get the segment with the given index.

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

    @property
    def prev(self):
        """Iteration: Previous iteration."""
        if self.number == 1:
            return None
        return self.dataset.iteration(self.number - 1)

    @property
    def next(self):
        """Iteration: Next iteration."""
        if self.number == self.dataset.num_iterations:
            return None
        return self.dataset.iteration(self.number + 1)

    def __eq__(self, other):
        return self.number == other.number and self.dataset is other.dataset

    def __repr__(self):
        return f'{self.__class__.__name__}({self.number}, {self.dataset})'


class WESTDataset:
    """A record of a WESTPA simulation.

    Parameters
    ----------
    h5filename : str, default 'west.h5'
        Pathname of a WESTPA HDF5 data file.
    segment_traj_loader : callable, optional
        A function to be used to load segment trajectories. The required
        signature and return type of `seg_traj_loader` are specified by the
        :class:`SegmentTrajectoryLoader` abstract base class.

    """

    def __init__(self, h5filename='west.h5', segment_traj_loader=None):
        self.h5file = WESTPAH5File(h5filename, 'r')
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
    def num_iterations(self):
        """int: Number of iterations in the simulation record."""
        return len(self.h5file['iterations'])

    @property
    def iterations(self):
        """Iterable[Iteration]: Sequence of iterations."""
        return (Iteration(n, self) for n in range(1, self.num_iterations + 1))

    def iteration(self, number):
        """Get the iteration with the given iteration number.

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

    def iteration_group(self, number):
        """Get the HDF5 group containing data for a given iteration.

        Parameters
        ----------
        number : int
            Iteration number (1-based).

        Returns
        -------
        h5py.Group
            The HDF5 group of the iteration indexed by `number`.

        """
        return self.h5file.get_iter_group(number)

    def __len__(self):
        return self.num_iterations

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.h5file.filename}')"


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


def trace(segment):
    """Return the trace (ancestral line) of a segment.

    Parameters
    ----------
    segment : Segment
        A segment of a :class:`WESTDataset`.

    Returns
    -------
    Iterable[Segment]
        The trace of `segment`.

    """
    if not isinstance(segment, Segment):
        raise TypeError('segment must be a Segment')

    segments = []
    while segment:
        segments.append(segment)
        segment = segment.parent

    return segments[::-1]
