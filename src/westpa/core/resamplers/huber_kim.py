import logging
import math

import numpy as np

from .base import ResamplerBase

logger = logging.getLogger(__name__)


class HuberKimResampler(ResamplerBase):
    """Implements the splitting and merging technique of Huber and Kim (1996). [1]_

    Parameters
    ----------
    adjust_counts : bool, default True
        Whether to adjust the number of walkers in occupied bins to exactly
        match the target count. This is a modification of the original
        Huber-Kim method, which only ensures that the number of walkers is
        close to the target count.
        Downward adjustments are made by iteratively merging the two
        lowest-weight walkers. Upward adjustments are made by iteratively
        splitting the highest-weight walker.
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    References
    ----------
    .. [1] G.A. Huber, S. Kim,
       Biophysical Journal, Volume 70, Issue 1, 1996, Pages 97-110, ISSN 0006-3495,
       https://doi.org/10.1016/S0006-3495(96)79552-8.

    """

    split_threshold = 2.0
    merge_cutoff = 1.0

    def __init__(self, adjust_counts=True, **kwargs):
        super().__init__(**kwargs)
        self.adjust_counts = adjust_counts

    def _split_by_weight(self, bin, ideal_weight):
        # Split walkers with weight > split_threshold * ideal_weight.
        index = bin.bisect_weights(self.split_threshold * ideal_weight)
        to_split = bin[index:]
        for segment in to_split:
            self.split_walker(bin, segment, m=math.ceil(segment.weight / ideal_weight))

    def _merge_by_weight(self, bin, ideal_weight):
        # Merge sets of walkers with combined weight <= merge_cutoff * ideal_weight.
        while True:
            cumulative_weight = bin.weights().cumsum()
            index = np.searchsorted(cumulative_weight, self.merge_cutoff * ideal_weight, side='right')
            to_merge = bin[:index]
            if len(to_merge) < 2:
                break
            self.merge_walkers(bin, to_merge, cumulative_weight[:index])

    def _adjust_count(self, bin, target_count):
        while len(bin) < target_count:
            logger.debug('adjusting counts by splitting')
            self.split_walker(bin, bin[-1])
        while len(bin) > target_count:
            logger.debug('adjusting counts by merging')
            self.merge_walkers(bin, bin[:2])

    def resample(self, bin, target_count):
        ideal_weight = bin.weight / target_count

        self._split_by_weight(bin, ideal_weight)
        self._merge_by_weight(bin, ideal_weight)

        if self.adjust_counts:
            self._adjust_count(bin, target_count)

        return bin
