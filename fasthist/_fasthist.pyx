# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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

from __future__ import division
import numpy, sys
cimport numpy, cython
from cpython.buffer cimport *
from cpython.mem cimport *
from numpy cimport * #PyArray_DATA, PyArray_TYPE

ctypedef fused real_numeric:
    numpy.int8_t
    numpy.int16_t
    numpy.int32_t
    numpy.int64_t
    numpy.uint8_t
    numpy.uint16_t
    numpy.uint32_t
    numpy.uint64_t    
    numpy.float32_t
    numpy.float64_t

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef histnd(values, binbounds, weights=1.0, out=None, binbound_check = True, ignore_out_of_range=False):
    '''Generate an N-dimensional PDF (or contribution to a PDF) from the given values.
    ``binbounds`` is a list of arrays of boundary values, with one entry for each
    dimension (``values`` must have as many columns as there are entries in ``binbounds``)
    ``weight``, if provided, specifies the weight each value contributes to the
    histogram; this may be a scalar (for equal weights for all values) or a vector of
    the same length as ``values`` (for unequal weights). If ``binbound_check`` is True, then
    the boundaries are checked for strict positive monotonicity; set to False to shave a few
    microseconds if you know your bin boundaries to be monotonically increasing. 
    '''
    
    if values.ndim != 2:
        values = numpy.atleast_2d(values)
        if values.ndim > 2:
            raise TypeError('values must be 2-D')
    
    cdef:
        Py_ssize_t npts = values.shape[0]
        Py_ssize_t ndim = values.shape[1] 
        int typecode = PyArray_TYPE(values) 
        
    if len(binbounds) != ndim:
        raise ValueError('number of sets of bin boundaries ({}) does not match dimensionality of data ({})'
                         .format(len(binbounds), values.shape[1]))

    if binbound_check:
        for idim in xrange(ndim):
            dq = numpy.diff(binbounds[idim])
            if (dq <= 0).any():
                raise ValueError('binbounds in dimension {} are not strictly monotonically increasing'.format(idim))

    # Prepare bin boundaries arrays
    _binbounds_vectors = numpy.empty((ndim,), numpy.object_)
    _nbounds   = numpy.empty((ndim,), numpy.uint32)
    for idim in range(ndim):
        _binbounds = numpy.require(binbounds[idim], values.dtype, 'C')
        _binbounds_vectors[idim] = _binbounds
        _nbounds[idim] = _binbounds.shape[0]

    # Prepare output array, if necessary        
    if out is None:
        _out = numpy.zeros([len(boundset)-1 for boundset in binbounds], numpy.float64)
    else:
        _out = out
        if _out.dtype != numpy.float64:
            raise TypeError('type of output array must be float64')
        if not _out.flags.writeable:
            raise TypeError('output is not writeable') 
        
    # Prepare weight array
    _weights = numpy.require(weights, numpy.float64, 'C')
    if _weights.shape == ():
        # scalar
        _weights = numpy.empty((len(values),), numpy.float64)
        _weights[:] = weights
    elif _weights.ndim > 1:
        raise TypeError('weight must be scalar or one dimensional')
    elif  _weights.shape[0] != values.shape[0]:
        raise TypeError('weights and values must be equal in length')
    
    # ugh
    if typecode == NPY_FLOAT32:      
        return _histnd[numpy.float32_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_FLOAT64:
        return _histnd[numpy.float64_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_INT8:
        return _histnd[numpy.int8_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_INT16:
        return _histnd[numpy.int16_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_INT32:
        return _histnd[numpy.int32_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_INT64:
        return _histnd[numpy.int64_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_UINT8:
        return _histnd[numpy.uint8_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_UINT16:
        return _histnd[numpy.uint16_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_UINT32:
        return _histnd[numpy.uint32_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    elif typecode == NPY_UINT64:
        return _histnd[numpy.uint64_t](values,
                                        _binbounds_vectors,
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        ignore_out_of_range,
                                        _out)
    else:
        raise TypeError('real floating-point or integer input required')
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef _histnd(real_numeric[:,:] values, object[:] binbounds, numpy.uint32_t* nbounds, double* weights,
             bint ignore_out_of_range, object output):
    '''Bin the values stored in the 2-D array ``values`` with corresponding weights ``weights``
    into the array of bins ``output``. The bin boundaries are specified in ``binbounds``, and the
    length of each set of boundaries is stored in ``nbounds``.
    
    Pre-conditions:
     * output is backed by 64-bit floating point storage
     * len(binbounds) == len(nbounds) == values.ndim
     '''
     
    cdef:
        Py_buffer outputview, bbview
        Py_ssize_t ndim = values.shape[1], npts = values.shape[0]
        Py_ssize_t idim, ipt, ibound
        char* outptr_bytes
        double* outptr
        real_numeric val, lb, ub
        real_numeric* boundbuf
        real_numeric** _binbounds
        bint store_value

    # Get pointers to the (contiguous) lists of bin boundaries in each dimension
    _binbounds = <real_numeric**> PyMem_Malloc(ndim*sizeof(real_numeric*))
    if not _binbounds:
        raise MemoryError()
    for idim in range(ndim):
        PyObject_GetBuffer(binbounds[idim], &bbview, PyBUF_SIMPLE)
        _binbounds[idim] = <real_numeric*> bbview.buf
        PyBuffer_Release(&bbview)

    # Get a view of our output array, so we can write directly into it 
    PyObject_GetBuffer(output, &outputview, PyBUF_STRIDED)

    try:
        with nogil:
            # loop over points
            for ipt in range(npts):
                outptr_bytes = <char*> outputview.buf
                store_value = True
                for idim in range(ndim):
                    val = values[ipt,idim]
                    for ibound in range(nbounds[idim]-1):
                        lb = _binbounds[idim][ibound]
                        ub = _binbounds[idim][ibound+1]
                        if val >= lb and val < ub:
                            outptr_bytes += outputview.strides[idim] * ibound
                            break
                    else:
                        if not ignore_out_of_range:
                            with gil:
                                raise ValueError('value {} at index {} out of bin boundaries in dimension {}'
                                                 .format(val,ipt,idim))
                        else:
                            store_value = False
                if store_value:
                    outptr = <double*> outptr_bytes
                    outptr[0] += weights[ipt]
        return output
    finally:
        PyMem_Free(_binbounds)
        PyBuffer_Release(&outputview)


