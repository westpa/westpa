import abc
import logging
import math
import secrets

import numpy as np

from westpa.core.we_driver import ConsistencyError

logger = logging.getLogger(__name__)


class ResamplerBase(abc.ABC):
    """Base class for resamplers.
    Subclasses must implement the :meth:`resample` method.

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

    _check = True

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
    def resample(self, bin_, target_count):
        """Resample the trajectories in a given bin.

        Parameters
        ----------
        bin_ : Bin
            Bin to resample.
        target_count : int
            Target number of trajectories for the bin.

        Returns
        -------
        bin_ : Bin
            Resampled bin.

        """
        ...

    def _split_by_threshold(self, bin):
        index = bin.bisect_weights(self.largest_allowed_weight, side='right')
        to_split = bin[index:]
        for segment in to_split:
            m = math.ceil(segment.weight / self.largest_allowed_weight)
            bin.split(segment, m=m)

    def _merge_by_threshold(self, bin):
        while True:
            index = bin.bisect_weights(self.smallest_allowed_weight)
            to_merge = bin[:index]
            if len(to_merge) < 2:
                return
            bin.merge(to_merge, rng=self.rng)

    def __call__(self, bin, target_count):
        if not bin:
            return bin  # skip empty bins

        bin_weight = bin.weight

        bin = self.resample(bin, target_count)

        if self._check:
            weights = bin.weights()
            if (weights <= 0).any():
                raise ConsistencyError('weights must be greater than 0')
            if not np.isclose(weights.sum(), bin_weight):  # TODO: Specify tolerance (atol, rtol).
                raise ConsistencyError('resampling must preserve the bin weight')

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
