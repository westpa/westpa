import numpy as np

from .base import ResamplerBase


def _resample_equal_weight(bin, counts, total_weight, target_count, new_weight=None):
    wtg_parent_ids = set.union(*(segment.wtg_parent_ids for segment in bin))
    new_weight = total_weight / target_count if new_weight is None else new_weight

    new_segments = set()
    for segment, count in zip(bin, counts):
        for _ in range(count):
            new_segment = segment.copy(weight=new_weight, wtg_parent_ids=wtg_parent_ids)
            new_segments.add(new_segment)

    bin.clear()
    bin |= new_segments

    return bin


class MultinomialResampler(ResamplerBase):
    """Implements the multinomial resampling technique.

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        weights = bin.weights()
        total_weight = weights.sum()

        counts = self.rng.multinomial(n=target_count, pvals=weights / total_weight)

        return _resample_equal_weight(bin, counts, total_weight, target_count)


class ResidualResampler(ResamplerBase):
    """Implements the residual resampling technique. (For details, see Algorithm 8.1
    of `this preprint <https://arxiv.org/abs/1806.00860>`_).

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        weights = bin.weights()
        total_weight = weights.sum()

        # See Algorithm 8.1 in https://arxiv.org/abs/1806.00860.
        nd = target_count * weights / total_weight
        nd_floor = np.floor(nd)
        delta = nd - nd_floor
        trials = round(delta.sum())
        if trials == 0:  # happens when bin.count == 1
            counts = [target_count]
        else:
            r = self.rng.multinomial(n=trials, pvals=delta / trials)
            counts = list(map(round, nd_floor + r))

        return _resample_equal_weight(bin, counts, total_weight, target_count)


class StratifiedResampler(ResamplerBase):
    """Implements the stratified resampling technique.
    (For details, see `this page <https://www.lancaster.ac.uk/stor-i-student-sites/martin-dimitrov/2021/05/14/resampling-techniques/>`_.)

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        cumulative_weight = bin.weights().cumsum()
        total_weight = cumulative_weight[-1]
        new_weight = total_weight / target_count

        random_values = self.rng.uniform(0, new_weight, size=target_count)
        random_values += np.linspace(0, total_weight - new_weight, num=target_count)

        indices = np.searchsorted(cumulative_weight, random_values)
        counts = np.bincount(indices)

        return _resample_equal_weight(bin, counts, total_weight, target_count, new_weight)


class SystematicResampler(ResamplerBase):
    """Implements the systematic resampling technique.
    (For details, see `this page <https://www.lancaster.ac.uk/stor-i-student-sites/martin-dimitrov/2021/05/14/resampling-techniques/>`_.)

    Parameters
    ----------
    **kwargs
        Keyword arguments to pass to the :class:`ResamplerBase` class.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def resample(self, bin, target_count):
        cumulative_weight = bin.weights().cumsum()
        total_weight = cumulative_weight[-1]
        new_weight = total_weight / target_count

        random_values = np.repeat(self.rng.uniform(0, new_weight), target_count)
        random_values += np.linspace(0, total_weight - new_weight, num=target_count)

        indices = np.searchsorted(cumulative_weight, random_values)
        counts = np.bincount(indices)

        return _resample_equal_weight(bin, counts, total_weight, target_count, new_weight)
