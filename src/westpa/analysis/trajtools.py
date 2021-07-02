import operator
import numpy as np

from concurrent.futures import ThreadPoolExecutor
from functools import cached_property, reduce
from tqdm import tqdm
from westpa.analysis.core import Walker, Trace

from typing import Callable, Iterable


class SegmentCollector:
    """An object that retrieves trajectory segments.

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

    def __init__(self, traj_descr, use_threads=False, max_workers=None,
                 show_progress=False):
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
        list
            The trajectory of each walker.

        """
        walkers = tuple(walkers)

        tqdm_kwargs = dict(
            desc='Retrieving segments',
            disable=(not self.show_progress),
            position=0,
            total=len(walkers),
        )

        if self.use_threads:
            with ThreadPoolExecutor(self.max_workers) as executor:
                segments = tqdm(executor.map(self.get_segment, walkers),
                                **tqdm_kwargs)
        else:
            segments = tqdm(map(self.get_segment, walkers), **tqdm_kwargs)

        return list(segments)


class Trajectory:
    """A data descriptor for walker and trace trajectories.

    Parameters
    ----------
    name : str
        Name of the :type:`Walker` and :type:`Trace` attribute to which to
        assign the descriptor.
    segment_getter : callable
        Function for getting a trajectory segment. Must take a single
        :type:`Walker` object as input and return a sequence representing
        the trajectory of the walker.
    cache_segments : bool, default True
        Whether to cache trajectory segments.

    See Also
    --------
    westpa.analysis.decorators.trajectory_segment
        Decorator that transforms a function for getting trajectory
        segments into a :type:`Trajectory` descriptor.

    """

    def __init__(self, name, segment_getter, cache_segments=True):
        for cls in Walker, Trace:
            if hasattr(cls, name):
                msg = f"class '{cls.__name__}' already has attribute '{name}'"
                raise AttributeError(msg)

        for cls in Walker, Trace:
            setattr(cls, name, self)

        self.name = name
        self.segment_getter = segment_getter
        self.cache_segments = cache_segments

    @cached_property
    def segment_collector(self):
        return SegmentCollector(self)

    @property
    def segment_getter(self):
        return self._segment_getter

    @segment_getter.setter
    def segment_getter(self, value):
        if not isinstance(value, Callable):
            raise TypeError('segment_getter must be callable')
        self._segment_getter = value

    @property
    def cache_segments(self):
        return self._cache_segments

    @cache_segments.setter
    def cache_segments(self, value):
        if not isinstance(value, bool):
            raise TypeError('cache_segments must be True or False')
        self._cache_segments = value

    @property
    def private_name(self):
        return '_' + self.name

    def __get__(self, instance, owner):
        if instance is None:
            return self

        if owner is Walker:
            if hasattr(instance, self.private_name):
                value = getattr(instance, self.private_name)
            else:
                value = self.segment_getter(instance)
                self.validate_segment(value)
                if self.cache_segments:
                    setattr(instance, self.private_name, value)
            return value

        if owner is Trace:
            segments = self.segment_collector.get_segments(instance)
            if isinstance(segments[0], np.ndarray):
                return np.concatenate(segments)
            return reduce(operator.concat, segments)

        msg = f'owner must be Walker or Trace, not {owner.__name__}'
        raise TypeError(msg)

    def __set__(self, instance, value):
        raise AttributeError("can't set attribute")

    def validate_segment(self, value):
        if not hasattr(value, '__getitem__'):
            msg = f"'{type(value).__name__}' object can't be concatenated"
            raise TypeError(msg)
