from __future__ import division
import numpy, sys
cimport numpy, cython
from cpython.buffer cimport *
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
cpdef histnd(values, binbounds, weights=1.0, out=None, binbound_check = True):
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
    _binbounds_ptrs = numpy.empty((ndim,), numpy.uintp)
    _nbounds   = numpy.empty((ndim,), numpy.uint32)
    for idim in range(ndim):
        _binbounds = numpy.require(binbounds[idim], values.dtype, 'C') 
        _binbounds_ptrs[idim] = <numpy.uintp_t> PyArray_DATA(_binbounds)
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
                                        <numpy.float32_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_FLOAT64:
        return _histnd[numpy.float64_t](values,
                                        <numpy.float64_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_INT8:
        return _histnd[numpy.int8_t](values,
                                        <numpy.int8_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_INT16:
        return _histnd[numpy.int16_t](values,
                                        <numpy.int16_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_INT32:
        return _histnd[numpy.int32_t](values,
                                        <numpy.int32_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_INT64:
        return _histnd[numpy.int64_t](values,
                                        <numpy.int64_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_UINT8:
        return _histnd[numpy.uint8_t](values,
                                        <numpy.uint8_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_UINT16:
        return _histnd[numpy.uint16_t](values,
                                        <numpy.uint16_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_UINT32:
        return _histnd[numpy.uint32_t](values,
                                        <numpy.uint32_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    elif typecode == NPY_UINT64:
        return _histnd[numpy.uint64_t](values,
                                        <numpy.uint64_t**> PyArray_DATA(_binbounds_ptrs),
                                        <numpy.uint32_t*> PyArray_DATA(_nbounds),
                                        <numpy.float64_t*> PyArray_DATA(_weights),
                                        _out)
    else:
        raise TypeError('real floating-point or integer input required')
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef _histnd(real_numeric[:,:] values, real_numeric** binbounds, numpy.uint32_t* nbounds, double* weights, object output):
    '''Bin the values stored in the 2-D array ``values`` with corresponding weights ``weights``
    into the array of bins ``output``. The bin boundaries are specified in ``binbounds``, and the
    length of each set of boundaries is stored in ``nbounds``.
    
    Pre-conditions:
     * output is backed by 64-bit floating point storage
     * len(binbounds) == len(nbounds) == values.ndim
     '''
     
    cdef:
        Py_buffer outputview
        Py_ssize_t ndim = values.shape[1], npts = values.shape[0]
        Py_ssize_t idim, ipt, ibound
        char* outptr_bytes
        double* outptr
        real_numeric val
        
    # Get a view of our output array, so we can write directly into it 
    PyObject_GetBuffer(output, &outputview, PyBUF_STRIDED)
    
    with nogil:
        # loop over points
        for ipt in range(npts):
            outptr_bytes = <char*> outputview.buf
            for idim in range(ndim):
                val = values[ipt,idim]
                for ibound in range(nbounds[idim]-1):
                    if binbounds[idim][ibound] <= val < binbounds[idim][ibound+1]:
                        outptr_bytes += outputview.strides[idim] * ibound
                        break
                else:
                    with gil:
                        raise ValueError('value {} at index {} out of bin boundaries in dimension {}'
                                         .format(val,ipt,idim))
            
            outptr = <double*> outptr_bytes
            outptr[0] += weights[ipt]
    return output

def _test_double(npts=1024*1024,loops=3):
    from time import time
    mine_times = [None]*loops
    theirs_times = [None]*loops
    
    # though 1.0 should be sufficient, weirdness in the boundary conditions
    # for numpy.digitize appears to introduce a discrepancy, so throw something
    # greater than 1.0 in for good measure
    binbounds = [[0,0.5,1,1.1] for x in xrange(3)]
    weights = numpy.random.rand(npts)
    #weights = numpy.ones((npts,), dtype=numpy.float64)
    for n in xrange(loops):
        testdat = numpy.random.rand(npts,3)
        mstart = time()
        mine = histnd(testdat, binbounds, weights=weights)
        mstop = time()
        tstart = time()
        theirs = numpy.histogramdd(testdat, binbounds, weights=weights)[0]
        tstop = time()
        mine_times[n] = mstop-mstart
        theirs_times[n] = tstop-tstart
        print mine
        print theirs
        errsum = numpy.abs(mine-theirs).sum()
        errsum_per_item = errsum / npts
        rel_err = errsum / numpy.abs(weights).sum()
        print 'sum of the absolute errors: {} ({} relative, {} per entry)'.format(errsum, rel_err, errsum_per_item)

    print 'mine, best of {}:   {}'.format(loops, min(mine_times))
    print 'theirs, best of {}: {}'.format(loops, min(theirs_times))
    
def _test_float(npts=1024*1024,loops=3):
    from time import time
    mine_times = [None]*loops
    theirs_times = [None]*loops
    
    binbounds = [[0,0.5,1,1.1] for x in xrange(3)]
    weights = numpy.random.rand(npts)
    #weights = numpy.ones((npts,), dtype=numpy.float64)
    for n in xrange(loops):
        testdat = numpy.require(numpy.random.rand(npts,3), numpy.float32)
        mstart = time()
        mine = histnd(testdat, binbounds, weights=weights)
        mstop = time()
        tstart = time()
        theirs = numpy.histogramdd(testdat, binbounds, weights=weights)[0]
        tstop = time()
        mine_times[n] = mstop-mstart
        theirs_times[n] = tstop-tstart
        print mine
        print theirs
        errsum = numpy.abs(mine-theirs).sum()
        errsum_per_item = errsum / npts
        rel_err = errsum / numpy.abs(weights).sum()
        print 'sum of the absolute errors: {} ({} relative, {} per entry)'.format(errsum, rel_err, errsum_per_item)

    print 'mine, best of {}:   {}'.format(loops, min(mine_times))
    print 'theirs, best of {}: {}'.format(loops, min(theirs_times))
  
def _test_uint(npts=1024*1024,loops=3):
    from time import time
    mine_times = [None]*loops
    theirs_times = [None]*loops
    
    binbounds = [[0,1,2,3,4] for x in xrange(3)]
    weights = numpy.random.rand(npts)
    #weights = numpy.ones((npts,), dtype=numpy.float64)
    for n in xrange(loops):
        testdat = numpy.require(numpy.random.randint(0,4,size=(npts,3)), numpy.uint16)
        mstart = time()
        mine = histnd(testdat, binbounds, weights=weights)
        mstop = time()
        tstart = time()
        theirs = numpy.histogramdd(testdat, binbounds, weights=weights)[0]
        tstop = time()
        mine_times[n] = mstop-mstart
        theirs_times[n] = tstop-tstart
        print mine
        print theirs
        errsum = numpy.abs(mine-theirs).sum()
        errsum_per_item = errsum / npts
        rel_err = errsum / numpy.abs(weights).sum()
        print 'sum of the absolute errors: {} ({} relative, {} per entry)'.format(errsum, rel_err, errsum_per_item)

    print 'mine, best of {}:   {}'.format(loops, min(mine_times))
    print 'theirs, best of {}: {}'.format(loops, min(theirs_times))
        
