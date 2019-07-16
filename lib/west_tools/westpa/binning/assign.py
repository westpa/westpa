# Copyright (C) 2013 Matthew C. Zwier, Joshua L. Adelman, and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

'''
Bin assignment for WEST simulations. This module defines "bin mappers" which take
vectors of coordinates (or rather, coordinate tuples), and assign each a definite
integer value identifying a bin. Critical portions are implemented in a Cython
extension module. 

A number of pre-defined bin mappers are available here:

  * :class:`RectilinearBinMapper`, for bins divided by N-dimensional grids
  * :class:`FuncBinMapper`, for functions which directly calculate bin assignments
    for a number of coordinate values. This is best used with C/Cython/Numba
    functions, or intellegently-tuned numpy-based Python functions.
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
(numpy.uint16) indicating a what bin that coordinate tuple falls into. The optional
``mask`` (a numpy bool array) specifies that some coordinates are to be skipped; this is used,
for instance, by the recursive (nested) bin mapper to minimize the number of calculations
required to definitively assign a coordinate tuple to a bin. Similarly, the optional
``output`` must be an integer (uint16) array of the same length as ``coords``, into which
assignments are written. The ``assign()`` function must return a reference to ``output``.
(This is used to avoid allocating many temporary output arrays in complex binning
scenarios.)

A user-defined bin mapper must also make an ``nbins`` property available, containing
the total number of bins within the mapper.

'''


import pickle as pickle
import hashlib, logging
import numpy

from . import _assign
from .bins import Bin

# All bin numbers are 16-bit unsigned ints, with one element (65525) reserved to
# indicate unknown or unassigned points. This allows up to 65,536 bins, making
# rate and flux matrices up to 32 GB (2**32 elements * 8 bytes). If you need more
# bins, change index_dtype here and index_dtype and index_t in _assign.pyx.
index_dtype = numpy.uint16
UNKNOWN_INDEX = 65535

# All coordinates are currently 32-bit floats. If you need 64-bit, change 
# coord_dtype here and coord_t in _assign.pyx.
coord_dtype = numpy.float32

from ._assign import output_map, apply_down, apply_down_argmin_across, rectilinear_assign 

log = logging.getLogger(__name__)

class BinMapper:
    hashfunc = hashlib.sha256

    def __init__(self):
        self.labels = None
        self.nbins = 0

    def construct_bins(self, type_=Bin):
        '''Construct and return an array of bins of type ``type``'''
        return numpy.array([type_() for _i in range(self.nbins)], dtype=numpy.object_)

    def pickle_and_hash(self):
        '''Pickle this mapper and calculate a hash of the result (thus identifying the
        contents of the pickled data), returning a tuple ``(pickled_data, hash)``.
        This will raise PickleError if this mapper cannot be pickled, in which case
        code that would otherwise rely on detecting a topology change must assume
        a topology change happened, even if one did not.
        '''

        pkldat = pickle.dumps(self, pickle.HIGHEST_PROTOCOL)
        hash = self.hashfunc(pkldat)
        return (pkldat, hash.hexdigest())

    def __repr__(self):
        return '<{} at 0x{:x} with {:d} bins>'.format(self.__class__.__name__, id(self), self.nbins or 0)

class NopMapper(BinMapper):
    '''Put everything into one bin.'''
    def __init__(self):
        super(NopMapper,self).__init__()
        self.nbins = 1
        self.labels = ['nop']

    def assign(self, coords, mask=None, output=None):
        if output is None:
            output = numpy.zeros((len(coords),), dtype=index_dtype)

        if mask is None:
            mask = numpy.ones((len(coords),), dtype=numpy.bool_)
        else:
            mask = numpy.require(mask, dtype=numpy.bool_)

        output[mask] = 0

class RectilinearBinMapper(BinMapper):
    '''Bin into a rectangular grid based on tuples of float values'''    
    def __init__(self, boundaries):
        super(RectilinearBinMapper,self).__init__()        
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
        self.labels = labels = []            
        for boundset in boundaries:
            boundarray = numpy.asarray(boundset,dtype=coord_dtype, order='C')
            db = numpy.diff(boundarray)
            if (db <= 0).any():
                raise ValueError('boundary set must be strictly monotonically increasing')
            self._boundaries.append(boundarray)
        self._boundlens = numpy.array([len(boundset) for boundset in self._boundaries], dtype=index_dtype)
        self.ndim = len(self._boundaries)
        self.nbins = numpy.multiply.accumulate([1] + [len(bounds)-1 for bounds in self._boundaries])[-1]

        _boundaries = self._boundaries
        binspace_shape = tuple(self._boundlens[:]-1)
        for index in numpy.ndindex(binspace_shape):
            bounds = [(_boundaries[idim][index[idim]], boundaries[idim][index[idim]+1]) for idim in range(len(_boundaries))]
            labels.append(repr(bounds))

    def assign(self, coords, mask=None, output=None):
        try:
            passed_coord_dtype = coords.dtype
        except AttributeError:
            coords = numpy.require(coords, dtype=coord_dtype)
        else:
            if passed_coord_dtype != coord_dtype:
                coords = numpy.require(coords, dtype=coord_dtype)

        if coords.ndim != 2:
            raise TypeError('coords must be 2-dimensional')

        if mask is None:
            mask = numpy.ones((len(coords),), dtype=numpy.bool_)
        elif len(mask) != len(coords):
            raise TypeError('mask [shape {}] has different length than coords [shape {}]'.format(mask.shape, coords.shape))

        if output is None:
            output = numpy.empty((len(coords),), dtype=index_dtype)
        elif len(output) != len(coords):
            raise TypeError('output has different length than coords')

        rectilinear_assign(coords, mask, output, self.boundaries, self._boundlens)

        return output

class PiecewiseBinMapper(BinMapper):
    '''Binning using a set of functions returing boolean values; if the Nth function
    returns True for a coordinate tuple, then that coordinate is in the Nth bin.'''

    def __init__(self, functions):
        self.functions = functions
        self.nbins = len(functions)
        self.index_dtype = numpy.min_scalar_type(self.nbins)
        self.labels = [repr(func) for func in functions]

    def assign(self, coords, mask=None, output=None):
        if output is None:
            output = numpy.zeros((len(coords),), dtype=index_dtype)

        if mask is None:
            mask = numpy.ones((len(coords),), dtype=numpy.bool_)
        else:
            mask = numpy.require(mask, dtype=numpy.bool_)

        coord_subset = coords[mask]
        fnvals = numpy.empty((len(coord_subset), len(self.functions)), dtype=index_dtype)
        for ifn, fn in enumerate(self.functions):
            rsl = numpy.apply_along_axis(fn,0,coord_subset)
            if rsl.ndim > 1:
                # this should work like a squeeze, unless the function returned something truly
                # stupid (e.g., a 3d array with at least two dimensions greater than 1), in which
                # case a broadcast error will occur
                fnvals[:,ifn] = rsl.flat
            else:
                fnvals[:,ifn] = rsl
        amask = numpy.require(fnvals.argmax(axis=1), dtype=index_dtype)
        output[mask] = amask
        return output    

class FuncBinMapper(BinMapper):
    '''Binning using a custom function which must iterate over input coordinate 
    sets itself.'''
    def __init__(self, func, nbins, args=None, kwargs=None):
        self.func = func
        self.nbins = nbins
        self.args = args or ()
        self.kwargs = kwargs or {}
        self.labels = ['{!r} bin {:d}'.format(func, ibin) for ibin in range(nbins)]

    def assign(self, coords, mask=None, output=None):
        try:
            passed_coord_dtype = coords.dtype
        except AttributeError:
            coords = numpy.require(coords, dtype=coord_dtype)
        else:
            if passed_coord_dtype != coord_dtype:
                coords = numpy.require(coords, dtype=coord_dtype)

        if coords.ndim != 2:
            raise TypeError('coords must be 2-dimensional')
        if mask is None:
            mask = numpy.ones((len(coords),), dtype=numpy.bool_)
        elif len(mask) != len(coords):
            raise TypeError('mask [shape {}] has different length than coords [shape {}]'.format(mask.shape, coords.shape))

        if output is None:
            output = numpy.empty((len(coords),), dtype=index_dtype)
        elif len(output) != len(coords):
            raise TypeError('output has different length than coords')

        self.func(coords, mask, output, *self.args, **self.kwargs)

        return output

class VectorizingFuncBinMapper(BinMapper):
    '''Binning using a custom function which is evaluated once for each (unmasked)
    coordinate tuple provided.'''
    def __init__(self, func, nbins, args=None, kwargs=None):
        self.func = func
        self.args = args or ()
        self.kwargs = kwargs or {}
        self.nbins = nbins
        self.index_dtype = numpy.min_scalar_type(self.nbins)
        self.labels = ['{!r} bin {:d}'.format(func, ibin) for ibin in range(nbins)]

    def assign(self, coords, mask=None, output=None):
        try:
            passed_coord_dtype = coords.dtype
        except AttributeError:
            coords = numpy.require(coords, dtype=coord_dtype)
        else:
            if passed_coord_dtype != coord_dtype:
                coords = numpy.require(coords, dtype=coord_dtype)

        if coords.ndim != 2:
            raise TypeError('coords must be 2-dimensional')
        if mask is None:
            mask = numpy.ones((len(coords),), dtype=numpy.bool_)
        elif len(mask) != len(coords):
            raise TypeError('mask [shape {}] has different length than coords [shape {}]'.format(mask.shape, coords.shape))

        if output is None:
            output = numpy.empty((len(coords),), dtype=index_dtype)
        elif len(output) != len(coords):
            raise TypeError('output has different length than coords')

        apply_down(self.func, self.args, self.kwargs, coords, mask, output)

        return output

class VoronoiBinMapper(BinMapper):
    '''A one-dimensional mapper which assigns a multidimensional pcoord to the
    closest center based on a distance metric. Both the list of centers and the
    distance function must be supplied.'''

    def __init__(self, dfunc, centers, dfargs=None, dfkwargs=None):
        self.dfunc = dfunc
        self.dfargs = dfargs or ()
        self.dfkwargs = dfkwargs or {}
        self.centers = numpy.asarray(centers)
        self.nbins = self.centers.shape[0]
        self.ndim = self.centers.shape[1]
        self.labels = ['center={!r}'.format(center) for center in self.centers]

        # Sanity check: does the distance map the centers to themselves?
        check = self.assign(self.centers)
        if (check != numpy.arange(len(self.centers))).any():
            raise TypeError('dfunc does not map centers to themselves')

    def assign(self, coords, mask=None, output=None):
        try:
            passed_coord_dtype = coords.dtype
        except AttributeError:
            coords = numpy.require(coords, dtype=coord_dtype)
        else:
            if passed_coord_dtype != coord_dtype:
                coords = numpy.require(coords, dtype=coord_dtype)

        if coords.ndim != 2:
            raise TypeError('coords must be 2-dimensional')
        if mask is None:
            mask = numpy.ones((len(coords),), dtype=numpy.bool_)
        elif len(mask) != len(coords):
            raise TypeError('mask [shape {}] has different length than coords [shape {}]'.format(mask.shape, coords.shape))

        if output is None:
            output = numpy.empty((len(coords),), dtype=index_dtype)
        elif len(output) != len(coords):
            raise TypeError('output has different length than coords')

        apply_down_argmin_across(self.dfunc, (self.centers,) + self.dfargs, self.dfkwargs, self.nbins,
                                 coords, mask, output)

        return output

class RecursiveBinMapper(BinMapper):
    '''Nest mappers one within another.'''

    def __init__(self, base_mapper, start_index=0):
        self.base_mapper = base_mapper
        self.nbins = base_mapper.nbins

        # Targets for recursion
        self._recursion_targets = {}

        # Which bins must we recurse into?
        self._recursion_map = numpy.zeros((self.base_mapper.nbins,), dtype=numpy.bool_)

        self.start_index = start_index

    @property
    def labels(self):
        for ilabel in range(self.base_mapper.nbins):
            if self._recursion_map[ilabel]:
                for label in self._recursion_targets[ilabel].labels:
                    yield label
            else:
                yield self.base_mapper.labels[ilabel]

    @property
    def start_index(self):
        return self._start_index

    @start_index.setter
    def start_index(self, new_index):
        self._start_index = new_index
        not_recursed = ~self._recursion_map
        n_not_recursed = not_recursed.sum()
        if n_not_recursed == self.nbins:
            self._output_map = numpy.arange(self._start_index, self._start_index + self.nbins, dtype=index_dtype)
        elif n_not_recursed > 0:
            # This looks like uninitialized access, but self._output_map is always set during __init__
            # (by self.start_index = 0, or whatever value was passed in), so this modifies the existing
            # set chosen above
            self._output_map[not_recursed] = numpy.arange(self._start_index, self._start_index + n_not_recursed,
                                                          dtype=index_dtype)
        else:
            # No un-replaced bins
            self._output_map = None

        n_own_bins = self.base_mapper.nbins - self._recursion_map.sum()
        startindex = self.start_index + n_own_bins
        for mapper in self._recursion_targets.values():
            mapper.start_index = startindex
            startindex += mapper.nbins

    def add_mapper(self, mapper, replaces_bin_at):
        '''Replace the bin containing the coordinate tuple ``replaces_bin_at`` with the
        specified ``mapper``.'''

        replaces_bin_at = numpy.require(replaces_bin_at, dtype=coord_dtype)
        if replaces_bin_at.ndim < 1:
            replaces_bin_at.shape = (1,1)
        elif replaces_bin_at.ndim < 2:
            replaces_bin_at.shape = (1,replaces_bin_at.shape[0])
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
        self.start_index = self.start_index        

    def assign(self, coords, mask=None, output=None):
        if mask is None:
            mask = numpy.ones((len(coords),), dtype=numpy.bool_)

        if output is None:
            output = numpy.empty((len(coords),), dtype=index_dtype)

        # mapping mask -- which output values come from our base
        # region set and therefore must be remapped            
        mmask = numpy.zeros((len(coords),), dtype=numpy.bool_)

        # Assign based on this mapper
        self.base_mapper.assign(coords, mask, output)

        # Which coordinates do we need to reassign, because they landed in
        # bins with embedded mappers?
        rmasks = {}
        for (rindex, mapper) in self._recursion_targets.items():
            omask = (output == rindex)
            mmask |= omask
            rmasks[rindex] = omask

        # remap output from our (base) mapper
        # omap may be None if every bin has a recursive mapper in it
        omap = self._output_map
        if omap is not None:
            output_map(output, omap, mask & ~mmask)

        # do any recursive assignments necessary
        for (rindex, mapper) in self._recursion_targets.items():
            mapper.assign(coords, mask&rmasks[rindex], output)

        return output
