from __future__ import print_function,division
import cython
import numpy
cimport numpy

# Coordinates are assumed to be 32-bit floats, and bin indices are
# defined to be 16-bit ints.

from numpy import uint16, float32
from numpy cimport uint16_t, float32_t

ctypedef numpy.float32_t _fptype  

ctypedef numpy.uint8_t bool_t
ctypedef float32_t coord_t    
ctypedef uint16_t index_t

bool_dtype = numpy.bool_
internal_bool_dtype = numpy.uint8
index_dtype = numpy.uint16
coord_dtype = numpy.float32

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef rectilinear_assign(numpy.ndarray[coord_t,ndim=2] coords,
                        numpy.ndarray[bool_t,ndim=1,cast=True] mask,
                        numpy.ndarray[index_t,ndim=1] output,
                        boundaries,
                        numpy.ndarray[index_t,ndim=1] boundlens):

    '''For bins delimited by sets boundaries on a rectilinear grid (``boundaries``),
    assign coordinates to bins, assuming C ordering of indices within the grid.
    ``boundlens`` is the number of boundaries in each dimension. 
    
    '''
    cdef:
        int icoord, idim, ibound, boundlen
        int ndim
        int ncoords = coords.shape[0]
        index_t index, stridefac
        coord_t bound
        coord_t cval

        numpy.ndarray[coord_t, ndim=1] boundvec
        numpy.ndarray[numpy.uintp_t, ndim=1] boundvecs
        coord_t* bvec
            
    # We assume greater locality across boundary vectors than across the
    # coordinate array. To exploit this locality, and only have to
    # relinquish/reacquire the GIL once, we extract pointers to the first
    # element of each boundary vector; we then release the GIL and go to
    # town on the entire data set.
    
    ndim = len(boundaries)
    boundvecs = numpy.empty((ndim,), dtype=numpy.uintp) 
    
    for 0 <= idim < ndim:
        boundvec = boundaries[idim]
        boundvecs[idim] = <numpy.uintp_t> &boundvec[0]
                
    with nogil:
        for icoord in range(ncoords):
            if not mask[icoord]:
                continue
            
            output[icoord] = 0
            stridefac = 1
    
            # backwards iteration needs signed values, so that the final != -1 works        
            for idim in range(ndim-1,-1,-1):
                found = 0
                cval = coords[icoord,idim]
                boundlen = boundlens[idim]
                bvec = <coord_t*> boundvecs[idim]
                
                if cval < bvec[0] or cval >= bvec[boundlen-1]:
                    with gil:
                        raise ValueError('coordinate value {} is out of bin space in dimension {}'.format(cval,idim))
                                    
                for ibound in range(1,boundlen):
                    if cval < bvec[ibound]:
                        index = ibound-1
                        break
                
                output[icoord] += index * stridefac
                stridefac *= boundlen-1
    return
        
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef testfunc(numpy.ndarray[coord_t, ndim=2] coords,
               numpy.ndarray[bool_t, ndim=1, cast=True]  mask,
               numpy.ndarray[index_t, ndim=1] output):
    cdef:
        index_t icoord
    
    for icoord in range(len(coords)):
        if mask[icoord]:
            if coords[icoord,0] < 0.5:
                output[icoord] = 0
            else:
                output[icoord] = 1
    return

# optimized function applications
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef apply_down(func,
                 args,
                 kwargs,
                 numpy.ndarray[coord_t, ndim=2] coords,
                 numpy.ndarray[bool_t, ndim=1, cast=True] mask,
                 numpy.ndarray[index_t, ndim=1] output):
    '''Apply func(coord, *args, **kwargs) to each input coordinate tuple,
    skipping any for which mask is false and writing results to output.'''
    cdef:
        Py_ssize_t i, n

    n = len(output)
    for i from 0 <= i < n:
        if mask[i]:
            output[i] = func(coords[i], *args, **kwargs)
    return

@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef apply_down_argmin_across(func,
                               args,
                               kwargs,
                               func_output_len,
                               numpy.ndarray[coord_t, ndim=2] coords,
                               numpy.ndarray[bool_t, ndim=1, cast=True] mask,
                               numpy.ndarray[index_t, ndim=1] output):
    '''Apply func(coord, *args, **kwargs) to each input coordinate tuple,
    skipping any for which mask is false and writing results to output.'''
    cdef:
        Py_ssize_t icoord, iout, ncoord, nout,
        coord_t _min
        index_t _argmin
        numpy.ndarray[coord_t, ndim=1] func_output
    
    nout = func_output_len
    func_output = numpy.empty((func_output_len,), dtype=coord_dtype)

    ncoord = len(coords)
    for icoord from 0 <= icoord < ncoord:
        if mask[icoord]:
            func_output = func(coords[icoord], *args, **kwargs)
            if len(func_output) != func_output_len:
                raise TypeError('function returned a vector of length {} (expected length {})'
                                .format(len(func_output), func_output_len))
            
            # find minimum value
            _min = func_output[0]
            _argmin = 0
            for iout from 1 <= iout < nout:
                if func_output[iout] < _min:
                    _min = func_output[iout]
                    _argmin = iout
                    
            output[icoord] = _argmin
    return
                
# optimized lookup table routine
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef output_map(numpy.ndarray[index_t, ndim=1] output,
                 numpy.ndarray[index_t, ndim=1] omap,
                 numpy.ndarray[bool_t, ndim=1, cast=True] mask):
    '''For each output for which mask is true, execute output[i] = omap[output[i]]'''

    cdef:
        Py_ssize_t i, n
        index_t o
        
    n = len(output)
    with nogil:
        for i from 0 <= i < n:
            if mask[i]:
                o = output[i]
                if o >= n:
                    with gil:
                        raise IndexError('value {} not available in output table'.format(o))
                output[i] = omap[o]
    return
