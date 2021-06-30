import functools
import operator
import numpy as np

from concurrent.futures import ThreadPoolExecutor
from westpa.analysis.core import Segment, Trace

from typing import Callable


class TrajectoryConcatenator:
    """An object that retrieves and concatenates segment trajectories.

    Parameters
    ----------
    segment_traj_getter : callable
        Function for getting a segment trajectory. Must take a single
        :type:`Segment` object as input and return a sequence.
    use_threads : bool, default False
        Whether to use a pool of threads to retrieve segment trajectories
        asynchronously. Setting this parameter to True may be may be
        useful when trajectory retrieval is an I/O bound task.
    max_workers : int, optional
        Maximum number of threads to use. The default value is specified
        in the :type:`ThreadPoolExecutor`
        `documentation <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor>`_.

    """

    def __init__(self, segment_traj_getter, use_threads=False,
                 max_workers=None):
        self._segment_traj_getter = segment_traj_getter
        self.use_threads = use_threads
        self.max_workers = max_workers

    @property
    def segment_traj_getter(self):
        return self._segment_traj_getter

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

    def join(self, segments):
        """Return the concatenated trajectory of a sequence of segments.

        Parameters
        ----------
        segments : Iterable[Segment]
            A sequence of segments (e.g., a trace).

        Returns
        -------
        sequence
            The concatenated trajectory of `segments`.

        """
        if self.use_threads:
            with ThreadPoolExecutor(self.max_workers) as executor:
                trajs = list(executor.map(self.segment_traj_getter, segments))
        else:
            trajs = list(map(self.segment_traj_getter, segments))
        if isinstance(trajs[0], np.ndarray):
            return np.concatenate(trajs)
        return functools.reduce(operator.concat, trajs)


class Trajectory:
    """A data descriptor for segment and trace trajectories.

    Parameters
    ----------
    name : str
        The name of the trajectory. When `cache_segments` is True, a
        segment trajectory will be cached as an attribute of the
        corresponding :type:`Segment` object, with attribute name
        ``'_' + name``. It is recommended that `name` be the name of the
        class attribute to which the descriptor is assigned.
    concatenator : TrajectoryConcatenator
        The object responsible for retrieving and concatenating segment
        trajectories.
    cache_segments : bool, default True
        Whether to cache segment trajectories.

    See Also
    --------
    westpa.analysis.decorators.segment_trajectory
        Decorator that transforms a function for getting segment
        trajectories into a :type:`Trajectory` descriptor.

    """

    def __init__(self, name, concatenator, cache_segments=True):
        self.name = name
        self.concatenator = concatenator
        self.cache_segments = cache_segments

    @property
    def private_name(self):
        return '_' + self.name

    def __get__(self, instance, owner):
        if instance is None:
            return self

        if owner is Segment:
            if hasattr(instance, self.private_name):
                value = getattr(instance, self.private_name)
            else:
                value = self.concatenator.segment_traj_getter(instance)
                self.validate(value)
                if self.cache_segments:
                    setattr(instance, self.private_name, value)
            return value

        if owner is Trace:
            return self.concatenator.join(instance)

        msg = f'owner must be Segment or Trace, not {owner.__name__}'
        raise TypeError(msg)

    def __set__(self, instance, value):
        raise AttributeError("can't set attribute")

    def validate(self, value):
        if not hasattr(value, '__getitem__'):
            msg = f"'{type(value).__name__}' object can't be concatenated"
            raise TypeError(msg)
