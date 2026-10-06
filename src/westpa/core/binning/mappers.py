from abc import ABC, abstractmethod
from functools import partial
from operator import attrgetter

import numpy as np
from sortedcontainers import SortedDict


class BinMapperBase(ABC):
    """Base class for bin mappers.
    Subclasses must implement the :meth:`assign` method and the :attr:`labels` property.

    Parameters
    ----------
    coord_getter : callable, optional
        Function that returns the coordinate time series to use for binning
        a segment. Must accept a segment and return a 2-D array. Defaults to
        ``attrgetter('pcoord')``.

    Attributes
    ----------
    labels : iterable of str
    nbins : int
        Number of bins mapped to.
    coord_getter : callable
        Function that returns the coordinate time series to use for binning
        a segment.

    Examples
    --------

    >>> import westpa
    >>> class TwoBinMapper(westpa.BinMapperBase):
    ...     @property
    ...     def labels(self):
    ...         return ['-', '+']
    ...     def assign(self, coords, weights, output):
    ...         output[:] = coords[:, 0] > 0
    ...         return output
    ...

    """

    def __init__(self, coord_getter=None):
        self.coord_getter = coord_getter or attrgetter('pcoord')

    @property
    @abstractmethod
    def labels(self):
        """Bin labels."""
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
    def assign(self, coords, weights, output):
        """Assign walkers to bins.

        Parameters
        ----------
        coords : numpy.ndarray of shape (n, d)
            Coordinates of each walker.
        weights : numpy.ndarray of shape (n,)
            Weight of each walker.
        output : numpy.ndarray of shape (n,)
            Array for storing output, initialized to -1.

        Returns
        -------
        output : numpy.ndarray of shape (n,)
            Bin assignment for each walker. Array elements must be in
            ``range(self.nbins)``.

        """
        ...

    def _get_coord(self, segment, coord_index):
        return self.coord_getter(segment)[coord_index]

    def __call__(self, segments, coord_index=-1):
        get_coord = partial(self._get_coord, coord_index=coord_index)

        coords = np.array(list(map(get_coord, segments)))
        weights = np.array(list(map(attrgetter('weight'), segments)))
        output = np.full(len(segments), -1)

        return self.assign(coords, weights, output)


class RecursiveBinMapper:
    """Nest bin mappers within one another.

    Parameters
    ----------
    base_mapper : BinMapper
        Bin mapper in which to nest other bin mappers.
    recursion_targets : Mapping[int, RecursiveBinMapper], optional
        Nested recursive mappers, keyed by `base_mapper` bin index.

    Attributes
    ----------
    nbins : int
        Total number of bins mapped to.
    labels : iterator of str
        Bin labels.
    base_mapper : BinMapper
        Base bin mapper.
    recursion_targets : Mapping[int, RecursiveBinMapper]
        Nested recursive mappers.

    Examples
    --------
    >>> import westpa
    >>> outer_mapper = westpa.RectilinearBinMapper([[0.0, 1.0, 2.0, 3.0]])
    >>> inner_mapper = westpa.RectilinearBinMapper([[1.0, 1.5, 2.0]])
    >>> rmapper = westpa.RecursiveBinMapper(outer_mapper)
    >>> rmapper.add_mapper(inner_mapper, replaces_bin=1)
    >>> rmapper.nbins
    4
    >>> list(rmapper.labels)
    ['[(0.0, 1.0)]', '[(2.0, 3.0)]', '[(1.0, 1.5)]', '[(1.5, 2.0)]']

    """

    def __init__(self, base_mapper, recursion_targets=None):
        self._base_mapper = base_mapper
        self._recursion_targets = SortedDict(recursion_targets or {})

    @property
    def base_mapper(self):
        return self._base_mapper

    @property
    def recursion_targets(self):
        return self._recursion_targets

    @property
    def labels(self):
        for idx, label in enumerate(self.base_mapper.labels):
            if idx not in self.recursion_targets:
                yield label
        for mapper in self.recursion_targets.values():
            yield from mapper.labels

    @property
    def nbins(self):
        return len(list(self.labels))

    def add_mapper(self, mapper, replaces_bin):
        """Replace the bin indexed by `replaces_bin` with the given bin mapper.

        Parameters
        ----------
        mapper : BinMapper
            Bin mapper with which to replace the bin indexed by `replaces_bin`.
        replaces_bin : int
            Index of the bin to replace.

        """
        not_recursed = [idx for idx in range(self.base_mapper.nbins) if idx not in self.recursion_targets]

        if replaces_bin < len(not_recursed):
            # replace a base mapper bin
            idx = not_recursed[replaces_bin]
            self.recursion_targets[idx] = RecursiveBinMapper(mapper)
        else:
            # replace a recursively embedded bin
            start_index = len(not_recursed)
            for nested_mapper in self.recursion_targets.values():
                if (idx := replaces_bin - start_index) in range(mapper.nbins):
                    nested_mapper.add_mapper(mapper, idx)
                    break
                start_index += nested_mapper.nbins

    def __call__(self, segments, coord_index=-1):
        assignments = np.full(len(segments), -1)
        base_assignments = self.base_mapper(segments, coord_index=coord_index)

        # non-recursive assignments
        output_map = {}
        for idx in range(self.base_mapper.nbins):
            if idx not in self.recursion_targets:
                output_map[idx] = len(output_map)
        mask = [idx in output_map for idx in base_assignments]
        assignments[mask] = [output_map[idx] for idx in base_assignments[mask]]

        # recursive assignments
        start_index = self.base_mapper.nbins - len(self.recursion_targets)
        for idx, mapper in self.recursion_targets.items():
            mask = base_assignments == idx
            if to_reassign := [segments[i] for i in np.flatnonzero(mask)]:
                assignments[mask] = mapper(to_reassign, coord_index=coord_index) + start_index
                start_index += mapper.nbins

        return assignments
