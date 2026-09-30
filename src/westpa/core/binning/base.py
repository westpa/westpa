from abc import ABC, abstractmethod
from operator import attrgetter

import numpy as np

from .bins import Bin


class BinMapperBase(ABC):
    """Base class for user-defined bin mappers.
    Subclasses must implement the :meth:`assign` method and the :attr:`labels` property.

    Parameters
    ----------
    coord_getter : callable, optional
        Function that returns a segment's binning coordinates. Must accept a
        segment and return a 2-D array (a time series). Defaults to
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

    UNKNOWN_INDEX = np.uint16(65535)  #: Indicates an unassigned segment.

    def __init__(self, coord_getter=None):
        self.coord_getter = coord_getter or attrgetter('pcoord')

    @property
    @abstractmethod
    def labels(self):
        pass

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

    def construct_bins(self):
        """Return a list of :attr:`nbins` empty bins.

        Returns
        -------
        bins : list of Bin
            List of empty bins, labeled according to :attr:`labels`.

        """
        return [Bin(label=label) for label in self.labels]

    @abstractmethod
    def assign(self, coords, weights, output):
        """Assign walkers to bins, given their coordinates and weights.

        Parameters
        ----------
        coords : 2-D numpy.ndarray
            Coordinates of each walker.
        weights : 1-D numpy.ndarray of dtype float64
            Weight of each walker.
        output : 1-D numpy.ndarray of dtype uint16
            Array for storing output, initialized to ``UNKNOWN_INDEX``.

        Returns
        -------
        output : 1-D numpy.ndarray of dtype uint16
            Bin assignment of each walker.

        """
        ...

    def _initial_coord(self, segment):
        return self.coord_getter(segment)[0]  # noqa

    def _final_coord(self, segment):
        return self.coord_getter(segment)[-1]  # noqa

    def __call__(self, segments, initial=False):
        get_coord = self._initial_coord if initial else self._final_coord

        coords = np.array(list(map(get_coord, segments)))
        weights = np.array(list(map(attrgetter('weight'), segments)))
        output = np.repeat(self.UNKNOWN_INDEX, len(segments))

        return self.assign(coords, weights, output)
