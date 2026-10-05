from collections import UserDict

import numpy as np


# Used to validate the 'nested_mappers' attribute.
# Ensures that keys are integers and 0 <= key < base_mapper.nbins.
class _NestedMappers(UserDict):

    def __init__(self, base_mapper_nbins, mapping=None):
        self.base_mapper_nbins = base_mapper_nbins
        super().__init__(mapping or {})

    def __setitem__(self, key, value):
        if isinstance(key, np.ndarray):
            key = key.item()
        if not isinstance(key, int):
            raise TypeError('keys must be integers, not ' + type(key).__name__)
        if key not in range(self.base_mapper_nbins):
            raise ValueError(f"keys must be in range({self.base_mapper_nbins})")
        super().__setitem__(key, value)


class RecursiveBinMapper:
    """Nest bin mappers within one another.

    Parameters
    ----------
    base_mapper : BinMapper
        Bin mapper in which to nest other bin mappers.
    nested_mappers : Mapping[int, BinMapper]
        Nested bin mappers, keyed by `base_mapper` bin index.

    Attributes
    ----------
    nbins : int
        Total number of bins mapped to.
    labels : iterator of str
        Bin labels.
    base_mapper : BinMapper
        Base bin mapper.
    nested_mappers : MutableMapping[int, BinMapper]
        Nested bin mappers.

    Examples
    --------
    >>> import westpa
    >>> bin_mapper = westpa.RecursiveBinMapper(
    ...     base_mapper=westpa.RectilinearBinMapper([[0.0, 1.0, 2.0, 3.0]]),
    ...     nested_mappers={
    ...         1: westpa.RectilinearBinMapper([[1.0, 1.25, 1.5, 1.75, 2.0]])
    ...     },
    ... )
    >>> bin_mapper.nbins
    6
    >>> list(bin_mapper.labels)
    ['[(0.0, 1.0)]',
     '[(2.0, 3.0)]',
     '[(1.0, 1.25)]',
     '[(1.25, 1.5)]',
     '[(1.5, 1.75)]',
     '[(1.75, 2.0)]']

    """

    def __init__(self, base_mapper, nested_mappers=None):
        self._base_mapper = base_mapper
        self._nested_mappers = _NestedMappers(base_mapper.nbins, nested_mappers)

    @property
    def base_mapper(self):
        return self._base_mapper

    @property
    def nested_mappers(self):
        return self._nested_mappers

    @property
    def labels(self):
        for idx, label in enumerate(self.base_mapper.labels):
            if idx not in self.nested_mappers:  # non-recursed
                yield label
        for mapper in self.nested_mappers.values():  # recursed
            yield from mapper.labels

    @property
    def nbins(self):
        return len(list(self.labels))

    def __call__(self, segments, coord_index=-1):
        assignments = np.full(len(segments), -1)
        base_assignments = self.base_mapper(segments, coord_index=coord_index)

        # non-recursed
        output_map = {}
        for idx in range(self.base_mapper.nbins):
            if idx not in self.nested_mappers:
                output_map[idx] = len(output_map)
        mask = [idx in output_map for idx in base_assignments]
        assignments[mask] = [output_map[idx] for idx in base_assignments[mask]]

        # recursed
        start_index = self.base_mapper.nbins - len(self.nested_mappers)
        for idx, mapper in self.nested_mappers.items():
            mask = base_assignments == idx
            if to_reassign := [segments[i] for i in np.flatnonzero(mask)]:
                assignments[mask] = mapper(to_reassign, coord_index=coord_index) + start_index
                start_index += mapper.nbins

        return assignments
