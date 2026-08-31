"""
Bin assignment for WEST simulations. This module defines "bin mappers" which take
vectors of coordinates (or rather, coordinate tuples), and assign each a definite
integer value identifying a bin. Critical portions are implemented in a Cython
extension module.

A number of pre-defined bin mappers are available here:

  * :class:`RectilinearBinMapper`, for bins divided by N-dimensional grids
  * :class:`FuncBinMapper`, for functions which directly calculate bin assignments
    for a number of coordinate values. This is best used with C/Cython/Numba
    functions, or intelligently tuned numpy-based Python functions.
  * :class:`VectorizingFuncBinMapper`, for functions which calculate a bin
    assignment for a single coordinate value. This is best used for arbitrary
    Python functions.
  * :class:`PiecewiseBinMapper`, for using a set of boolean-valued functions, one
    per bin, to determine assignments. This is likely to be much slower than a
    `FuncBinMapper` or `VectorizingFuncBinMapper` equipped with an appropriate
    function, and its use is discouraged.

One "super-mapper" is available, for assembling more complex bin spaces from
simpler components:

  * :class:`RecursiveBinMapper`, for nesting one set of bins within another.

Users are also free to implement their own mappers. A bin mapper must implement, at
least, an ``assign(coords, mask=None, output=None)`` method, which is responsible
for mapping each of the vector of coordinate tuples ``coords`` to an integer
(np.uint16) indicating a what bin that coordinate tuple falls into. The optional
``mask`` (a numpy bool array) specifies that some coordinates are to be skipped; this is used,
for instance, by the recursive (nested) bin mapper to minimize the number of calculations
required to definitively assign a coordinate tuple to a bin. Similarly, the optional
``output`` must be an integer (uint16) array of the same length as ``coords``, into which
assignments are written. The ``assign()`` function must return a reference to ``output``.
(This is used to avoid allocating many temporary output arrays in complex binning
scenarios.)

A user-defined bin mapper must also make an ``nbins`` property available, containing
the total number of bins within the mapper.

"""

import hashlib
import logging
import pickle

import numpy as np

from .bins import Bin
from ._assign import output_map, apply_down, apply_down_argmin_across, rectilinear_assign

# All bin numbers are 16-bit unsigned ints, with one element (65525) reserved to
# indicate unknown or unassigned points. This allows up to 65,536 bins, making
# rate and flux matrices up to 32 GB (2**32 elements * 8 bytes). If you need more
# bins, change index_dtype here and index_dtype and index_t in _assign.pyx.
index_dtype = np.uint16
UNKNOWN_INDEX = 65535

# All coordinates are currently 32-bit floats. If you need 64-bit, change
# coord_dtype here and coord_t in _assign.pyx.
coord_dtype = np.float32

log = logging.getLogger(__name__)


class BinMapper:
    """Base class for bin mappers. Subclasses must implement the :meth:`map`
    method, as well as provide values for the :attr:`nbins` and :attr:`labels`
    attributes.

    Attributes
    ----------
    nbins : int
        Total number of bins to which the mapper assigns trajectories.
    labels : list of str
        Label for each bin.

    Methods
    -------
    map

    """

    def __init__(self):
        self.labels = None
        self.nbins = 0

    def construct_bins(self):
        """Return a list of empty bins."""
        if self.labels:
            return [Bin(label=label) for label in self.labels]
        else:
            return [Bin() for _ in range(self.nbins)]

    def pickle_and_hash(self):
        """Pickle this mapper and calculate a hash of the result (thus identifying the
        contents of the pickled data), returning a tuple ``(pickled_data, hash)``.
        This will raise PickleError if this mapper cannot be pickled, in which case
        code that would otherwise rely on detecting a topology change must assume
        a topology change happened, even if one did not.
        """
        pkldat = pickle.dumps(self, pickle.HIGHEST_PROTOCOL)
        hash = hashlib.sha256(pkldat)
        return (pkldat, hash.hexdigest())

    def __repr__(self):
        return '<{} at 0x{:x} with {:d} bins>'.format(self.__class__.__name__, id(self), self.nbins or 0)

    def _assign(self, coords, mask, output):
        raise NotImplementedError()

    def assign(self, coords, mask=None, output=None):
        # Validate arguments and pass them to _assign(coords, mask, output),
        # which handles the actual bin assignment.
        try:
            passed_coord_dtype = coords.dtype
        except AttributeError:
            coords = np.require(coords, dtype=coord_dtype)
        else:
            if passed_coord_dtype != coord_dtype:
                coords = np.require(coords, dtype=coord_dtype)

        if coords.ndim != 2:
            raise TypeError('coords must be 2-dimensional')

        if mask is None:
            mask = np.ones((len(coords),), dtype=np.bool_)
        elif len(mask) != len(coords):
            raise TypeError('mask [shape {}] has different length than coords [shape {}]'.format(mask.shape, coords.shape))

        if output is None:
            output = np.empty((len(coords),), dtype=index_dtype)
        elif len(output) != len(coords):
            raise TypeError('output has different length than coords')

        self._assign(coords, mask, output)

        return output

    def __call__(self, segments):
        coords = np.array(list(map(lambda seg: seg.pcoord[-1], segments)))
        assignments = self.assign(coords)

        bins = self.construct_bins()
        for segment, idx in zip(segments, assignments):
            bins[idx].add(segment)

        return bins


class NopMapper(BinMapper):
    """Put everything into one bin."""

    def __init__(self):
        super().__init__()
        self.nbins = 1
        self.labels = ['nop']

    def _assign(self, coords, mask, output):
        output[mask] = 0


class RectilinearBinMapper(BinMapper):
    """Bin into a rectangular grid.

    Parameters
    ----------
    boundaries : iterable of array_like
        Bin boundaries along each progress coordinate dimension.

    """

    def __init__(self, boundaries):
        super().__init__()
        self._boundaries = None
        self._boundlens = None
        self.ndim = 0
        self.nbins = 0

        # the setter function below handles all of the required wrangling
        self.boundaries = boundaries

    @property
    def boundaries(self):
        return self._boundaries

    @boundaries.setter
    def boundaries(self, boundaries):
        del self._boundaries, self.labels
        self._boundaries = []
        self.labels = []
        for boundset in boundaries:
            boundarray = np.asarray(boundset, dtype=coord_dtype)
            if not boundarray.ndim == 1:
                raise ValueError("items in 'boundaries' must be 1-D arrays")
            db = np.diff(boundarray)
            if (db <= 0).any():
                raise ValueError('bin boundaries must be monotonically increasing along each dimension')
            self._boundaries.append(boundarray)
        self._boundlens = np.array([len(boundset) for boundset in self._boundaries], dtype=index_dtype)
        self.ndim = len(self._boundaries)
        self.nbins = np.multiply.accumulate([1] + [len(bounds) - 1 for bounds in self._boundaries])[-1]

        _boundaries = self._boundaries
        binspace_shape = tuple(self._boundlens[:] - 1)
        for index in np.ndindex(binspace_shape):
            label = (
                '['
                + ', '.join(
                    f'({boundarray[index[idim]]!s}, {boundarray[index[idim] + 1]!s})' for idim, boundarray in enumerate(_boundaries)
                )
                + ']'
            )
            self.labels.append(label)

    def _assign(self, coords, mask, output):
        rectilinear_assign(coords, mask, output, self.boundaries, self._boundlens)


class PiecewiseBinMapper(BinMapper):
    """Binning using a set of functions returning boolean values; if the Nth function
    returns True for a coordinate tuple, then that coordinate is in the Nth bin."""

    def __init__(self, functions):
        self.functions = functions
        self.nbins = len(functions)
        self.index_dtype = np.min_scalar_type(self.nbins)
        self.labels = [str(func) for func in functions]

    def _assign(self, coords, mask, output):
        coord_subset = coords[mask]
        fnvals = np.empty((len(coord_subset), len(self.functions)), dtype=index_dtype)
        for ifn, fn in enumerate(self.functions):
            rsl = np.apply_along_axis(fn, 0, coord_subset)
            if rsl.ndim > 1:
                # this should work like a squeeze, unless the function returned something truly
                # stupid (e.g., a 3d array with at least two dimensions greater than 1), in which
                # case a broadcast error will occur
                fnvals[:, ifn] = rsl.flat
            else:
                fnvals[:, ifn] = rsl
        amask = np.require(fnvals.argmax(axis=1), dtype=index_dtype)
        output[mask] = amask


class FuncBinMapper(BinMapper):
    """Binning using a custom function which must iterate over input coordinate
    sets itself."""

    def __init__(self, func, nbins, args=None, kwargs=None):
        self.func = func
        self.nbins = nbins
        self.args = args or ()
        self.kwargs = kwargs or {}
        self.labels = ['{!r} bin {:d}'.format(func, ibin) for ibin in range(nbins)]

    def _assign(self, coords, mask, output):
        self.func(coords, mask, output, *self.args, **self.kwargs)


class VectorizingFuncBinMapper(BinMapper):
    """Binning using a custom function which is evaluated once for each (unmasked)
    coordinate tuple provided."""

    def __init__(self, func, nbins, args=None, kwargs=None):
        self.func = func
        self.args = args or ()
        self.kwargs = kwargs or {}
        self.nbins = nbins
        self.index_dtype = np.min_scalar_type(self.nbins)
        self.labels = ['{!r} bin {:d}'.format(func, ibin) for ibin in range(nbins)]

    def _assign(self, coords, mask, output):
        apply_down(self.func, self.args, self.kwargs, coords, mask, output)


class VoronoiBinMapper(BinMapper):
    """Assign progress coordinate points to the closest center based on a
    distance metric. Both the list of centers and the distance function must
    be supplied.

    Parameters
    ----------
    dfunc : callable
        Distance function. It must accept arguments ``(x, ys)`` and return
        a 1-D array containing the distance of each point in ``ys`` to the
        point ``x``.
    centers : 2-D array_like
        Voronoi sites.
    dfargs : tuple, optional
        Optional arguments to pass to `dfunc`.
    dfkwargs : Mapping[str, Any], optional
        Optional keyword arguments to pass to `dfunc`.

    """

    def __init__(self, dfunc, centers, dfargs=None, dfkwargs=None):
        self.dfunc = dfunc
        self.dfargs = dfargs or ()
        self.dfkwargs = dfkwargs or {}
        self.centers = np.asarray(centers)
        self.nbins = self.centers.shape[0]
        self.ndim = self.centers.shape[1]
        self.labels = ['center={!r}'.format(center) for center in self.centers]

        # Sanity check: does the distance map the centers to themselves?
        check = self.assign(self.centers)
        if (check != np.arange(len(self.centers))).any():
            raise TypeError('dfunc does not map centers to themselves')

    def _assign(self, coords, mask, output):
        apply_down_argmin_across(self.dfunc, (self.centers,) + self.dfargs, self.dfkwargs, self.nbins, coords, mask, output)


class RecursiveBinMapper(BinMapper):
    """Nest mappers one within another.

    Parameters
    ----------
    base_mapper : BinMapper
        Base mapper within which to nest other bin mappers.
    start_index : int, default 0
        Initial bin index.

    Methods
    -------
    add_mapper

    """

    def __init__(self, base_mapper, start_index=0):
        self.base_mapper = base_mapper
        self.nbins = base_mapper.nbins

        # Targets for recursion
        self._recursion_targets = {}

        # Which bins must we recurse into?
        self._recursion_map = np.zeros((self.base_mapper.nbins,), dtype=np.bool_)

        self.start_index = start_index

    @property
    def labels(self):
        for ilabel in np.flatnonzero(~self._recursion_map):
            yield self.base_mapper.labels[ilabel]
        for mapper in self._recursion_targets.values():
            for label in mapper.labels:
                yield label

    @property
    def start_index(self):
        return self._start_index

    @start_index.setter
    def start_index(self, new_index):
        self._start_index = new_index
        not_recursed = ~self._recursion_map
        n_not_recursed = not_recursed.sum()
        if n_not_recursed == self.nbins:
            self._output_map = np.arange(self._start_index, self._start_index + self.nbins, dtype=index_dtype)
        elif n_not_recursed > 0:
            # This looks like uninitialized access, but self._output_map is always set during __init__
            # (by self.start_index = 0, or whatever value was passed in), so this modifies the existing
            # set chosen above
            self._output_map[not_recursed] = np.arange(self._start_index, self._start_index + n_not_recursed, dtype=index_dtype)
        else:
            # No un-replaced bins
            self._output_map = None

        n_own_bins = self.base_mapper.nbins - self._recursion_map.sum()
        startindex = self.start_index + n_own_bins
        for mapper in self._recursion_targets.values():
            mapper.start_index = startindex
            startindex += mapper.nbins

    def add_mapper(self, mapper, replaces_bin_at):
        """Replace the bin containing the coordinate tuple `replaces_bin_at` with the
        specified `mapper`.

        Parameters
        ----------
        mapper : BinMapper
            Bin mapper with which to replace the bin containing
            `replaces_bin_at`.
        replaces_bin_at : array_like
            Coordinate tuple indicating the bin to replace.

        """
        replaces_bin_at = np.require(replaces_bin_at, dtype=coord_dtype)
        if replaces_bin_at.ndim < 2:
            replaces_bin_at = np.atleast_2d(replaces_bin_at)
        elif replaces_bin_at.ndim > 2 or replaces_bin_at.shape[1] > 1:
            raise TypeError('a single coordinate vector is required')

        self.nbins += mapper.nbins - 1

        ibin = self.base_mapper.assign(replaces_bin_at)[0]
        log.debug('replacing bin {!r} containing {!r} with {!r}'.format(ibin, replaces_bin_at, mapper))
        if self._recursion_map[ibin]:
            # recursively add; this doesn't change anything for us except our
            # total bin count, which has been accounted for above
            self._recursion_targets[ibin].add_mapper(mapper, replaces_bin_at[0])
        else:
            # replace a bin on our mapper
            self._recursion_map[ibin] = True
            mapper = RecursiveBinMapper(mapper)
            self._recursion_targets[ibin] = mapper

        # we have updated our list of recursed bins, so set our own start index to trigger a recursive
        # reassignment of mappers' output values
        # Note that we're reordering the recursion targets based on outer bin numbers (dict keys) first,
        # so the order the mappers were added no longer matters...
        self._recursion_targets = {k: self._recursion_targets[k] for k in sorted(self._recursion_targets)}
        self.start_index = self.start_index

    def _assign(self, coords, mask, output):
        # mapping mask -- which output values come from our base
        # region set and therefore must be remapped
        mmask = np.zeros((len(coords),), dtype=np.bool_)

        # Assign based on this mapper
        self.base_mapper.assign(coords, mask, output)

        # Which coordinates do we need to reassign, because they landed in
        # bins with embedded mappers?
        rmasks = {}
        for rindex, mapper in self._recursion_targets.items():
            omask = output == rindex
            mmask |= omask
            rmasks[rindex] = omask

        # remap output from our (base) mapper
        # omap may be None if every bin has a recursive mapper in it
        omap = self._output_map
        if omap is not None:
            output_map(output, omap, mask & ~mmask)

        # do any recursive assignments necessary
        for rindex, mapper in self._recursion_targets.items():
            mapper.assign(coords, mask & rmasks[rindex], output)
