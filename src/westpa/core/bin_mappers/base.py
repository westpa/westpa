from abc import ABC, abstractmethod
from functools import partial
from operator import attrgetter

import numpy as np


class BinMapperBase(ABC):
    """Base class for bin mappers.
    Subclasses must implement the :meth:`map` method and the :attr:`labels` property.

    Parameters
    ----------
    coord_getter : callable, optional
        Function that returns the coordinate time series to use for binning a
        segment. Must accept a segment and return a 2-D array. Defaults to
        ``attrgetter('pcoord')``.

    Attributes
    ----------
    labels : iterable of str
        Bin labels.
    nbins : int
        Number of bins mapped to.
    coord_getter : callable
        Function that returns a segment's binning coordinates.

    Examples
    --------

    >>> import westpa
    >>> class TwoBinMapper(westpa.BinMapperBase):
    ...     @property
    ...     def labels(self):
    ...         return ['-', '+']
    ...     def map(self, coords, weights, output):
    ...         output[:] = coords[:, 0] > 0
    ...         return output
    ...

    """

    def __init__(self, coord_getter=None):
        self.coord_getter = coord_getter or attrgetter('pcoord')

    @property
    @abstractmethod
    def labels(self):  # noqa (black/flake8 one liner impasse)
        ...

    @property
    def nbins(self):
        return len(list(self.labels))

    @property
    def coord_getter(self):
        return self._coord_getter

    @coord_getter.setter
    def coord_getter(self, value):
        if not callable(value):
            raise TypeError("'coord_getter' must be callable")
        self._coord_getter = value

    @abstractmethod
    def map(self, coords, weights, output):
        """Map walkers to bins, given their coordinates and weights.

        Parameters
        ----------
        coords : 2-D numpy.ndarray
            Coordinates of each walker.
        weights : 1-D numpy.ndarray of float
            Weight of each walker.
        output : 1-D numpy.ndarray of int
            Array for storing output, initialized to -1.

        Returns
        -------
        output : 1-D numpy.ndarray of float
            Bin assignment of each walker.

        """
        ...

    def _get_coord(self, segment, coord_index):
        return self.coord_getter(segment)[coord_index]

    def __call__(self, segments, coord_index=-1):
        get_coord = partial(self._get_coord, coord_index=coord_index)

        coords = np.array(list(map(get_coord, segments)))
        weights = np.array(list(map(attrgetter('weight'), segments)))
        output = np.full(len(segments), -1)

        return self.map(coords, weights, output)
