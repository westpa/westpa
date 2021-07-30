import concurrent.futures as cf
import operator
import numpy as np

# Required for fast concatenation of MDTraj trajectories.
try:
    import mdtraj
except ImportError:
    mdtraj = None

from functools import cached_property, partial, reduce
from tqdm import tqdm
from westpa.analysis.core import Walker, Trace

from typing import Callable


class SegmentCollector:
    """An object that manages the retrieval of trajectory segments.

    Parameters
    ----------
    traj_descr : Trajectory
        A trajectory descriptor.
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

    def __init__(self, traj_descr, use_threads=False, max_workers=None, show_progress=False):
        self.traj_descr = traj_descr
        self.use_threads = use_threads
        self.max_workers = max_workers
        self.show_progress = show_progress

    @property
    def traj_descr(self):
        return self._traj_descr

    @traj_descr.setter
    def traj_descr(self, value):
        if not isinstance(value, Trajectory):
            msg = f'traj_descr must be an instance of {Trajectory}'
            raise TypeError(msg)
        self._traj_descr = value

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

    def get_segment(self, walker):
        """Return the trajectory of a given walker.

        Parameters
        ----------
        walker : Walker
            A walker.

        Returns
        -------
        sequence
            The trajectory of `walker`.

        """
        return self.traj_descr.__get__(walker, Walker)

    def get_segments(self, walkers):
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

        if self.use_threads:
            with cf.ThreadPoolExecutor(self.max_workers) as executor:
                future_to_key = {executor.submit(self.get_segment, walker): key for key, walker in enumerate(walkers)}
                futures = list(tqdm(cf.as_completed(future_to_key), **tqdm_kwargs))
                futures.sort(key=future_to_key.get)
                segments = (future.result() for future in futures)
        else:
            segments = tqdm(map(self.get_segment, walkers), **tqdm_kwargs)

        return list(segments)


class Trajectory:
    """A data descriptor for walker and trace trajectories.

    Parameters
    ----------
    fget : callable
        Function for getting a trajectory segment. Must take a single
        :type:`Walker` object as input and return a sequence representing
        the trajectory of the walker.
    name : str, optional
        Name of the :type:`Walker` and :type:`Trace` attribute to which
        to assign the trajectory descriptor. If not provided, `name` will
        default to the function name of `fget`.
    concatenator : callable, optional
        Function for concatenating trajectories. Must take a sequence of
        trajectories as input and return their concatenation. The default
        `concatenator` is :func:`concatenate`.
    cache_segments : bool, default True
        Whether to cache trajectory segments.

    See Also
    --------
    :func:`trajectory_segment`
        Decorator that transforms a function for getting trajectory
        segments into a :type:`Trajectory` descriptor.

    """

    def __init__(self, fget=None, *, name=None, concatenator=None, cache_segments=True):
        if fget is None:
            return partial(self.__init__, name=name, concatenator=concatenator, cache_segments=cache_segments)

        if name is None:
            name = fget.__name__

        self.fget = fget
        self.name = name
        self.concatenator = concatenator
        self.cache_segments = cache_segments

        # Attach self to Walker and Trace classes.
        for cls in Walker, Trace:
            if hasattr(cls, name):
                msg = f"class '{cls.__name__}' already has attribute '{name}'"
                raise AttributeError(msg)
        for cls in Walker, Trace:
            setattr(cls, name, self)

    @property
    def private_name(self):
        """str: Name of the :type:`Walker` instance attribute used for
        caching segments.

        """
        return '_' + self.name

    @cached_property
    def segment_collector(self):
        """SegmentCollector: Segment retrieval manager."""
        return SegmentCollector(self)

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
    def cache_segments(self):
        """bool: Whether to cache trajectory segments."""
        return self._cache_segments

    @cache_segments.setter
    def cache_segments(self, value):
        if not isinstance(value, bool):
            raise TypeError('cache_segments must be True or False')
        self._cache_segments = value

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

    def __get__(self, instance, owner):
        if instance is None:
            return self

        if owner is Walker:
            if hasattr(instance, self.private_name):
                value = getattr(instance, self.private_name)
            else:
                value = self.fget(instance)
                self._validate_segment(value)
                if self.cache_segments:
                    setattr(instance, self.private_name, value)
            return value

        if owner is Trace:
            segments = self.segment_collector.get_segments(instance)
            return self.concatenator(segments)

        msg = f'owner must be Walker or Trace, not {owner.__name__}'
        raise TypeError(msg)

    def __set__(self, instance, value):
        raise AttributeError("can't set attribute")

    def __call__(self, arg):
        if isinstance(arg, Walker):
            return self.__get__(arg, Walker)
        if isinstance(arg, Trace):
            return self.__get__(arg, Trace)
        raise TypeError('argument must be a Walker or Trace')

    def _validate_segment(self, value):
        if not hasattr(value, '__getitem__'):
            msg = f"'{type(value).__name__}' object can't be concatenated"
            raise TypeError(msg)


def trajectory_segment(fget=None, *, cache=True):
    """Transform a function for getting a trajectory segment into a
    trajectory attribute of the same name.

    Parameters
    ----------
    fget : callable
        Function for getting a trajectory segment. Must take a single
        :type:`Walker` object as input and return a sequence.
    cache : bool, default True
        Whether to cache trajectory segments.

    Returns
    -------
    Trajectory
        The newly created trajectory attribute. Equivalent to
        ``getattr(Walker, fget.__name__)`` and
        ``getattr(Trace, fget.__name__)``.

    """
    return Trajectory(fget, cache_segments=cache)


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
