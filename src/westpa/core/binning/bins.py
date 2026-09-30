import logging
import operator
from collections.abc import MutableSet, Sequence
from functools import reduce

import numpy as np
from sortedcontainers import SortedSet

logger = logging.getLogger(__name__)

EPS = np.finfo(np.float64).eps


class Bin(MutableSet, Sequence):
    """Mutable set of segments, sorted in increasing order of weight.
    The order is automatically maintained as the bin is updated.

    Parameters
    ----------
    segments : iterable of Segment, optional
        Initial set of segments.
    label : str, optional
        Descriptive label for the bin.

    Attributes
    ----------
    label : str
        Bin label.
    weight : float
        Total weight of all the segments in the bin.

    Examples
    --------

    Create a bin containing three segments:

    >>> import westpa
    >>> segments = [westpa.Segment(weight=weight) for weight in [0.3, 0.1, 0.2]]
    >>> bin_ = westpa.Bin(segments)
    >>> bin_.weights()
    array([0.1, 0.2, 0.3])

    Get the segment with the lowest or highest weight:

    >>> bin_[0]
    <Segment n_iter=None, seg_id=None, weight=0.1, parent_id=None at 0x1056d2990>
    >>> bin_[-1]
    <Segment n_iter=None, seg_id=None, weight=0.3, parent_id=None at 0x16a4f5580>

    """

    def __init__(self, segments=None, label=None):
        self._segments = SortedSet(segments, key=operator.attrgetter('weight'))
        self._label = str(label) if label else ''

    def __repr__(self):
        return f'<{type(self).__name__} label={self.label!r}, count={len(self)}, weight={self.weight} at {hex(id(self))}>'

    # The next five methods are required to implement MutableSet.
    def __contains__(self, segment):
        """Return True if the segment is in the bin, False otherwise."""
        return segment in self._segments

    def __iter__(self):
        """Return an iterator over the segments in the bin."""
        return iter(self._segments)

    def __len__(self):
        """Return the number of segments in the bin."""
        return len(self._segments)

    def add(self, segment):
        """Add a segment to the bin."""
        self._segments.add(segment)

    def discard(self, segment):
        """Remove a segment from the bin. Do not raise an exception if absent."""
        self._segments.discard(segment)

    # MutableSet doesn't provide update() or difference_update(), used by we_driver.
    def update(self, *others):
        self._segments.update(*others)

    def difference_update(self, *others):
        self._segments.difference_update(*others)

    # Default clear() mixin method is slow (calls pop() repeatedly).
    def clear(self):
        """Empty the bin."""
        self._segments.clear()

    # Implement Sequence.
    def __getitem__(self, index):
        """Return the segment at the given index. Supports slicing."""
        return self._segments[index]

    @property
    def label(self):
        return self._label

    @property
    def weight(self):
        return sum(map(self._segments.key, self))

    def weights(self):
        """Return the segment weights in increasing order.

        Returns
        -------
        weights : 1-D numpy.ndarray
            Sorted array of segment weights.

        """
        return np.array(list(map(self._segments.key, self)))

    def reweight(self, new_weight):
        """Reweight the bin by scaling the segment weights.

        Parameters
        ----------
        new_weight : float
            New bin weight after reweighting. Must be between 0 and 1.

        Returns
        -------
        self : Bin
            Reweighted bin.

        """
        if not (0 <= new_weight <= 1):
            raise ValueError("'new_weight' must be between 0 and 1")

        if len(self) == 0:
            if new_weight > 0:
                raise ValueError('cannot reweight empty bin')
            return self

        ratio = new_weight / self.weight
        segments = [segment.copy(weight=ratio * segment.weight) for segment in self]
        self._segments = SortedSet(segments, key=operator.attrgetter('weight'))

        return self

    def bisect_weights(self, w, side='left'):
        """Find the index where `w` should be inserted in ``self.weights()`` to maintain sorted order.

        Parameters
        ----------
        w : float
            Value to insert.
        side : {'left', 'right'}, optional
            If 'left', return the insertion point before (to the left of) any
            existing entries of `w`. If 'right', return the insertion point
            after (to the right of) any existing entries of `w`.

        Returns
        -------
        index : int
            Insertion point for `w`.

        """
        match side:
            case 'left':
                return self._segments.bisect_key_left(w)
            case 'right':
                return self._segments.bisect_key_right(w)
            case _:
                raise ValueError("'side' must be either 'left' or 'right'")

    def split(self, segment, m=2):
        """Split a segment into two or more copies.

        Parameters
        ----------
        segment : Segment
            Segment to split. This segment is removed from the bin and
            replaced by the returned segments.
        m : int, default 2
            Number of copies to split `segment` into.

        Returns
        -------
        new_segments : list of Segment
            `m` copies of `segment`, each with weight ``segment.weight / m``.

        """
        if not isinstance(m, int):
            raise TypeError("'m' must be an integer")
        if not m >= 2:
            raise ValueError("'m' must be greater than or equal to 2")

        new_weight = segment.weight / m
        new_segments = [segment.copy(weight=new_weight) for _ in range(m)]

        self.remove(segment)
        self.update(new_segments)

        return new_segments

    def merge(self, segments, cumulative_weight=None, rng=None):
        """Merge multiple segments into a copy of a single segment. The surviving
        (copied) segment is chosen randomly according to weight.

        Parameters
        ----------
        segments : iterable of Segment
            Segments to merge. These segments are removed from the bin and
            replaced by the returned segment.
        cumulative_weight : 1-D array_like, optional
            Cumulative sum of the segment weights. If not passed, the value
            is computed by this method.
        rng : numpy.random.Generator, int, or sequence of int, optional
            Pseudorandom number generator (PRNG) to use, or a seed to
            initialize the PRNG. Integer values must be nonnegative.
            Defaults to ``numpy.random.default_rng()``.

        Returns
        -------
        new_segment : Segment
            Copy of the surviving segment, with weight equal to the combined
            weight of `segments`.

        """
        rng = np.random.default_rng(rng)

        segments = list(segments)
        weights = np.array(list(map(operator.attrgetter('weight'), segments)))

        if cumulative_weight is None:
            cumulative_weight = weights.cumsum()

        idx = np.digitize(rng.uniform(0, cumulative_weight[-1]), cumulative_weight)

        new_segment = segments[idx].copy(
            weight=cumulative_weight[-1],
            wtg_parent_ids=reduce(operator.or_, (segment.wtg_parent_ids for segment in segments)),
        )

        self.difference_update(segments)
        self.add(new_segment)

        return new_segment
