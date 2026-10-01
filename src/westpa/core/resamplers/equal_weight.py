import operator
from functools import reduce

import numpy as np

from .base import ResamplerBase


def _resample(bin, counts, new_weight):
    wtg_parent_ids = reduce(operator.or_, (segment.wtg_parent_ids for segment in bin))

    new_segments = set()
    for segment, count in zip(bin, counts):
        for _ in range(count):
            new_segment = segment.copy(weight=new_weight, wtg_parent_ids=wtg_parent_ids)
            new_segments.add(new_segment)

    bin.clear()
    bin.update(new_segments)

    return bin


class MultinomialResampler(ResamplerBase):
    """Multinomial resampling.

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    _check = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        weights = bin.weights()
        total_weight = weights.sum()

        counts = self.rng.multinomial(n=target_count, pvals=weights / total_weight)

        return _resample(bin, counts, total_weight / target_count)


class ResidualResampler(ResamplerBase):
    """Residual resampling.

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    _check = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        weights = bin.weights()
        total_weight = weights.sum()

        # See Algorithm 8.1 of https://arxiv.org/abs/1806.00860.
        nd = target_count * weights / total_weight
        nd_floor = np.floor(nd)
        delta = nd - nd_floor
        trials = round(delta.sum())
        if trials == 0:  # happens when len(bin) == 1
            counts = [target_count]
        else:
            r = self.rng.multinomial(n=trials, pvals=delta / trials)
            counts = list(map(round, nd_floor + r))

        return _resample(bin, counts, total_weight / target_count)


class StratifiedResampler(ResamplerBase):
    """Stratified resampling.

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    _check = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        cumulative_weight = bin.weights().cumsum()
        total_weight = cumulative_weight[-1]
        new_weight = total_weight / target_count

        offsets = np.linspace(0, total_weight - new_weight, num=target_count)
        random_vals = self.rng.uniform(0, new_weight, size=target_count) + offsets

        counts = np.bincount(np.searchsorted(cumulative_weight, random_vals))

        return _resample(bin, counts, new_weight)


class SystematicResampler(ResamplerBase):
    """Systematic resampling.

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    _check = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        cumulative_weight = bin.weights().cumsum()
        total_weight = cumulative_weight[-1]
        new_weight = total_weight / target_count

        offsets = np.linspace(0, total_weight - new_weight, num=target_count)
        random_vals = self.rng.uniform(0, new_weight) + offsets

        counts = np.bincount(np.searchsorted(cumulative_weight, random_vals))

        return _resample(bin, counts, new_weight)
