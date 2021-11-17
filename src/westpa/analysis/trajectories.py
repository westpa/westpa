import concurrent.futures as cf
import functools
import operator
import numpy as np

# Required for fast concatenation of MDTraj trajectories.
try:
    import mdtraj
except ImportError:
    mdtraj = None

from functools import partial, reduce
from tqdm import tqdm
from westpa.analysis.core import Walker, Trace

from typing import Callable


class Trajectory:
    """A callable that returns the trajectory of a walker or trace.

    Parameters
    ----------
    fget : callable
        Function for retrieving a single trajectory segment. Must take a
        :type:`Walker` object as its first argument and return a sequence
        (e.g., an MDTraj Trajectory or a NumPy array) representing the
        trajectory of the walker.
    concatenator : callable, optional
        Function for concatenating trajectories. Must take a sequence of
        trajectories as input and return their concatenation. The default
        `concatenator` is :func:`concatenate`.

    """

    def __init__(self, fget=None, *, concatenator=None):
        if fget is None:
            return partial(self.__init__, concatenator=concatenator)

        self.fget = fget
        self.concatenator = concatenator

        self._segment_collector = SegmentCollector(self)

    @property
    def segment_collector(self):
        """SegmentCollector: Segment retrieval manager."""
        return self._segment_collector

    @property
    def fget(self):
        """callable: Function for getting trajectory segments."""
        return self._fget

    @fget.setter
    def fget(self, value):
        if not isinstance(value, Callable):
            raise TypeError('fget must be callable')
        self._fget = value

    @property
    def concatenator(self):
        """callable: Function for concatenating trajectories."""
        return self._concatenator

    @concatenator.setter
    def concatenator(self, value):
        if value is None:
            value = concatenate
        elif not isinstance(value, Callable):
            raise TypeError('concatenator must be callable')
        self._concatenator = value

    def __call__(self, obj, **kwargs):
        if isinstance(obj, Walker):
            value = self.fget(obj, **kwargs)
            self._validate_segment(value)
            return value
        if isinstance(obj, Trace):
            segments = self.segment_collector.get_segments(obj.walkers, **kwargs)
            return self.concatenator(segments)
        raise TypeError('argument must be a Walker or Trace')

    def _validate_segment(self, value):
        if not hasattr(value, '__getitem__'):
            msg = f"'{type(value).__name__}' object can't be concatenated"
            raise TypeError(msg)


class SegmentCollector:
    """An object that manages the retrieval of trajectory segments.

    Parameters
    ----------
    trajectory : Trajectory
        The trajectory to which the segment collector is attached.
    use_threads : bool, default False
        Whether to use a pool of threads to retrieve trajectory segments
        asynchronously. Setting this parameter to True may be may be
        useful when segment retrieval is an I/O bound task.
    max_workers : int, optional
        Maximum number of threads to use. The default value is specified
        in the :type:`ThreadPoolExecutor`
        `documentation <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor>`_.
    show_progress : bool, default False
        Whether to show a progress bar when retrieving multiple segments.

    """

    def __init__(self, trajectory, use_threads=False, max_workers=None, show_progress=False):
        self.trajectory = trajectory
        self.use_threads = use_threads
        self.max_workers = max_workers
        self.show_progress = show_progress

    @property
    def trajectory(self):
        return self._trajectory

    @trajectory.setter
    def trajectory(self, value):
        if not isinstance(value, Trajectory):
            msg = f'trajectory must be an instance of {Trajectory}'
            raise TypeError(msg)
        self._trajectory = value

    @property
    def use_threads(self):
        return self._use_threads

    @use_threads.setter
    def use_threads(self, value):
        if not isinstance(value, bool):
            raise TypeError('use_threads must be True or False')
        self._use_threads = value

    @property
    def max_workers(self):
        return self._max_workers

    @max_workers.setter
    def max_workers(self, value):
        if value is None:
            self._max_workers = None
            return
        if value <= 0:
            raise ValueError('max_workers must be greater than 0')
        self._max_workers = value

    @property
    def show_progress(self):
        return self._show_progress

    @show_progress.setter
    def show_progress(self, value):
        if not isinstance(value, bool):
            raise ValueError('show_progress must be True or False')
        self._show_progress = value

    def get_segments(self, walkers, **kwargs):
        """Return the trajectories of a group of walkers.

        Parameters
        ----------
        walkers : Iterable[Walker]
            A group of walkers.

        Returns
        -------
        list of sequences
            The trajectory of each walker.

        """
        walkers = tuple(walkers)

        tqdm_kwargs = dict(desc='Retrieving segments', disable=(not self.show_progress), position=0, total=len(walkers),)

        get_segment = functools.partial(self.trajectory, **kwargs)
        if self.use_threads:
            with cf.ThreadPoolExecutor(self.max_workers) as executor:
                future_to_key = {executor.submit(get_segment, walker): key for key, walker in enumerate(walkers)}
                futures = list(tqdm(cf.as_completed(future_to_key), **tqdm_kwargs))
                futures.sort(key=future_to_key.get)
                segments = (future.result() for future in futures)
        else:
            segments = tqdm(map(get_segment, walkers), **tqdm_kwargs)

        return list(segments)


def concatenate(trajectories):
    """Return the concatenation of a sequence of trajectories.

    Parameters
    ----------
    trajectories : sequence of sequences
        A sequence of trajectories.

    Returns
    -------
    sequence
        The concatenation of `trajectories`.

    """
    if isinstance(trajectories[0], np.ndarray):
        return np.concatenate(trajectories)
    if mdtraj and isinstance(trajectories[0], mdtraj.Trajectory):
        return trajectories[0].join(trajectories[1:], check_topology=False)
    return reduce(operator.concat, trajectories)
