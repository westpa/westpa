import abc
import logging
import math
import operator
import secrets

import numpy as np

from westpa.core.we_driver import ConsistencyError

logger = logging.getLogger(__name__)


class ResamplerBase(abc.ABC):
    """Base class for resamplers. Subclasses must implement the :meth:`resample` method.

    Parameters
    ----------
    rng : numpy.random.Generator, int, or sequence of int, optional
        Psuedorandom number generator (PRNG) to use, or a seed for
        initializing the PRNG. Integer values must be nonnegative.
        Defaults to ``numpy.random.default_rng()``.
    smallest_allowed_weight : float, default 1e-310
        Minimum weight threshold.
    largest_allowed_weight : float, default 1.0
        Maximum weight threshold.
    thresholds : bool, default True
        Whether to enforce the `smallest_allowed_weight` and
        `largest_allowed_weight` thresholds. If True, this is done after
        calling :meth:`resample`.

    Attributes
    ----------
    rng : numpy.random.Generator
        Psuedorandom number generator.
    smallest_allowed_weight : float
        Minimum weight threshold.
    largest_allowed_weight : float
        Maximum weight threshold.
    thresholds : bool
        Whether to enforce the weight thresholds.

    """

    def __init__(
        self,
        rng=None,
        smallest_allowed_weight=1e-310,
        largest_allowed_weight=1.0,
        thresholds=True,
    ):
        if rng is None:
            seed = secrets.randbits(128)
            rng = np.random.default_rng(seed)
            logger.info(f'rng=default_rng({seed=})')
        else:
            rng = np.random.default_rng(rng)

        if not (0 < smallest_allowed_weight < 1):
            raise ValueError("'smallest_allowed_weight' must be between 0 and 1")
        if not (smallest_allowed_weight < largest_allowed_weight <= 1):
            raise ValueError("'largest_allowed_weight' must be between 'smallest_allowed_weight' and 1")

        self.rng = rng
        self.smallest_allowed_weight = smallest_allowed_weight
        self.largest_allowed_weight = largest_allowed_weight
        self.thresholds = thresholds

    @abc.abstractmethod
    def resample(self, bin, target_count):
        """Resample the walkers in a given bin.

        Parameters
        ----------
        bin : Bin
            Bin to resample.
        target_count : int
            Target number of walkers for the bin.

        Returns
        -------
        bin : Bin
            Resampled bin.

        """
        ...

    @staticmethod
    def split_walker(bin, segment, m=2):
        """Split a walker into two or more copies.

        This method modifies `bin` by replacing the input `segment` with the
        output `new_segments`.

        Parameters
        ----------
        bin : Bin
            Bin containing the walker.
        segment : Segment
            Walker to split.
        m : int, default 2
            Number of copies to split `segment` into.

        Returns
        -------
        new_segments : set of Segment
            New segments created by splitting `segment`.

        """
        if not isinstance(m, int):
            raise TypeError("'m' must be an integer")
        if not m >= 2:
            raise ValueError("'m' must be greater than or equal to 2")

        new_weight = segment.weight / m
        new_segments = {segment.copy(weight=new_weight) for _ in range(m)}

        bin.remove(segment)
        bin |= new_segments

        return new_segments

    def merge_walkers(self, bin, segments, cumulative_weight=None):
        """Merge multiple walkers into a single walker. The surviving walker
        is chosen randomly according to weight.

        This method modifies `bin` by replacing the input `segments` with the
        output `new_segment`.

        Parameters
        ----------
        bin : Bin
            Bin containing the walkers.
        segments : iterable of Segment
            Walkers to merge.
        cumulative_weight : 1-D array-like, optional
            Cumulative sum of the walker weights. If not passed, the value will be
            computed by this function.

        Returns
        -------
        new_segment : Segment
            New segment created by merging `segments`.

        """
        segments = list(segments)
        weights = np.array(list(map(operator.attrgetter('weight'), segments)))

        if cumulative_weight is None:
            cumulative_weight = weights.cumsum()

        idx = np.digitize(self.rng.uniform(0, cumulative_weight[-1]), cumulative_weight)

        new_segment = segments[idx].copy(
            weight=cumulative_weight[-1],
            wtg_parent_ids=set.union(*(segment.wtg_parent_ids for segment in segments)),
        )

        bin -= segments
        bin.add(new_segment)

        return new_segment

    def _split_by_threshold(self, bin):
        index = bin.bisect_weights(self.largest_allowed_weight, side='right')
        to_split = bin[index:]
        for segment in to_split:
            m = math.ceil(segment.weight / self.largest_allowed_weight)
            self.split_walker(bin, segment, m=m)

    def _merge_by_threshold(self, bin):
        while True:
            index = bin.bisect_weights(self.smallest_allowed_weight)
            to_merge = bin[:index]
            if len(to_merge) < 2:
                return
            self.merge_walkers(bin, to_merge)

    def __call__(self, bin, target_count):
        if not bin:
            return bin

        initial_weight = bin.weight

        bin = self.resample(bin, target_count)

        weights = bin.weights()
        if (weights <= 0).any():
            raise ConsistencyError('weights must be greater than 0')
        if not math.isclose(weights.sum(), initial_weight, abs_tol=1e-12):  # TODO: What should this tolerance be?
            raise ConsistencyError('resampling must preserve the total weight of the bin')

        if self.thresholds:
            self._split_by_threshold(bin)
            self._merge_by_threshold(bin)
            for segment in bin:
                if not (self.smallest_allowed_weight <= segment.weight <= self.largest_allowed_weight):
                    logger.warning(
                        'Unable to fulfill weight threshold conditions for %s. The given threshold range is likely too small.',
                        segment,
                    )

        return bin
