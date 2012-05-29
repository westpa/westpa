from __future__ import division
import numpy, sys
cimport numpy, cython 

@cython.boundscheck(False)
@cython.wraparound(False)
def hist(values, binbounds, out = None,
         double cweight=1.0,
         binbound_check = True
         ):
    '''
    Generate a PDF (or a contribution to a PDF) from the given values.
    '''
    
    if binbound_check:
        tdx = numpy.diff(binbounds)
        if (tdx <= 0).any():
            raise ValueError('binbounds is not strictly monotonically increasing')
    
    cdef numpy.uint64_t ival, nval
    cdef numpy.int32_t ibb, nbb, ibb0, ibbu, ibbl
    cdef double val, minbound, maxbound, pivot_up, pivot_down
    
    cdef numpy.ndarray[numpy.float64_t,ndim=1] _values = numpy.array(values, numpy.float64, copy=False)
    cdef numpy.ndarray[numpy.float64_t,ndim=1] _binbounds = numpy.array(binbounds, numpy.float64, copy=False)
    cdef numpy.ndarray[numpy.float64_t,ndim=1] _histout

    nval = _values.shape[0]
    nbb = _binbounds.shape[0]
    minbound = _binbounds[0]
    maxbound = _binbounds[_binbounds.shape[0]-1]
            
    if out is None:
        _histout = numpy.zeros((nbb-1,), numpy.float64)
    else:
        if binbounds.shape[0] != out.shape[0] + 1:
            raise TypeError('incompatible dimensions for binbounds and out')        
        _histout = numpy.array(out, dtype=numpy.float64, copy=False)
    
    with nogil:
        ibb = nbb >> 1
        for ival in xrange(nval):
            val = _values[ival]
            if val < minbound or val > maxbound:
                with gil:
                    raise ValueError('element {:d} outside bin range'.format(ival))
            
            if val >= _binbounds[ibb] and val < _binbounds[ibb+1]:
                # found in current bin
                pass
            else:
                for ibb in xrange(0,nbb-1):                
                    if val >= _binbounds[ibb] and val < _binbounds[ibb+1]:
                        break

#                ibb0 = ibb
#                pivot_up = _binbounds[ibb0+1]
#                pivot_down = _binbounds[ibb]
#                
#                if val >= pivot_up:
#                    # need to search higher
#                    while ibb < nbb-2:
#                        if val >= _binbounds[ibb] and val < _binbounds[ibb]:
#                            break
#                        ibb+=1
#                    else:
#                        # val == maxbound                        
#                        if val > maxbound:
#                            with gil:
#                                raise ValueError('element {:d} outside bin range'.format(ival))
#                elif val < pivot_down:
#                    # need to search lower
#                    while ibb >= 1:
#                        if val >= _binbounds[ibb] and val < _binbounds[ibb+1]:
#                            break
#                        ibb -= 1
#                    else:
#                        if val < minbound:
#                            with gil:
#                                raise ValueError('element {:d} outside bin range'.format(ival))
#                else:
#                    with gil:
#                        raise AssertionError('impossible branch')
                
                
            
            _histout[ibb] += cweight #* (_binbounds[ibb+1] - _binbounds[ibb])
                                
    return _histout

@cython.boundscheck(False)
@cython.wraparound(False)
def hist2d(values, binbounds, out = None,
           double cweight=1.0,
           binbound_check = True
           ):
    '''
    Generate a PDF (or a contribution to a PDF) from the given values.
    '''
    
    if len(binbounds) != 2:
        raise ValueError('binbounds must be a list of [xbinbounds, ybinbounds]')
    else:
        binbounds_x = binbounds[0]
        binbounds_y = binbounds[1]
    
    if binbound_check:
        dx = numpy.diff(binbounds_x)
        if (dx <= 0).any():
            raise ValueError('binbounds in x are not strictly monotonically increasing')
        
        dy = numpy.diff(binbounds_y)
        if (dy <= 0).any():
            raise ValueError('binbounds in y are not strictly monotonically increasing')

    values = numpy.array(values, dtype=numpy.float64, copy=False)
    if values.ndim != 2 or 2 not in values.shape:
        raise TypeError('values must be convertible to an (N,2) array')
    elif values.shape[1] != 2:
        values.transpose()
            
    cdef numpy.uint64_t ival, nval
    cdef numpy.int32_t ibbx, ibb0x, nbbx, ibboffx
    cdef numpy.int32_t ibby, ibb0y, nbby, ibboffy
    cdef double val_x, val_y
    cdef double minbound_x, maxbound_x, minbound_y, maxbound_y
    
    cdef numpy.ndarray[numpy.float64_t,ndim=2] _values = values
    cdef numpy.ndarray[numpy.float64_t,ndim=1] _binbounds_x = binbounds_x
    cdef numpy.ndarray[numpy.float64_t,ndim=1] _binbounds_y = binbounds_y    
    cdef numpy.ndarray[numpy.float64_t,ndim=2] _histout

    nval = _values.shape[0]
    nbbx = _binbounds_x.shape[0]
    nbby = _binbounds_y.shape[0]
    minbound_x = _binbounds_x[0]
    maxbound_x = _binbounds_x[nbbx-1]
    minbound_y = _binbounds_y[0]
    maxbound_y = _binbounds_y[nbby-1]
            
    if out is None:
        _histout = numpy.zeros((nbbx-1, nbby-1), numpy.float64)
    elif out.ndim != 2:
        raise TypeError('output array must be of dimension 2')
    else:
        if (nbbx, nbby) != (out.shape[0] + 1, out.shape[1] + 1): 
            raise TypeError('incompatible dimensions for binbounds and out')        
        _histout = numpy.array(out, dtype=numpy.float64, copy=False)
    
    with nogil:
        for ival in xrange(nval):
            val_x = _values[ival,0]
            val_y = _values[ival,1]
            
            if val_x < minbound_x or val_x > maxbound_x or val_y < minbound_y or val_y > maxbound_y:
                with gil:
                    raise ValueError('element {:d} outside bin range'.format(ival))
            
            for ibbx in xrange(0,nbbx-1):                
                if val_x >= _binbounds_x[ibbx] and val_x < _binbounds_x[ibbx+1]:
                    break
            
            for ibby in xrange(0,nbby-1):
                if val_y >= _binbounds_y[ibby] and val_y < _binbounds_y[ibby+1]:
                    break
                                            
            _histout[ibbx,ibby] += cweight * ( (_binbounds_x[ibbx+1] - _binbounds_x[ibbx]) 
                                             * (_binbounds_y[ibby+1] - _binbounds_y[ibby]) )
                
    return _histout
