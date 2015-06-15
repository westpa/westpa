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

from __future__ import print_function,division
import cython
import numpy
import warnings
cimport numpy

ctypedef numpy.uint16_t index_t
ctypedef numpy.float64_t weight_t
ctypedef numpy.uint8_t bool_t
ctypedef numpy.int64_t seg_id_t
ctypedef numpy.uint_t uint_t # 32 bits on 32-bit systems, 64 bits on 64-bit systems

cdef double NAN = numpy.nan 

weight_dtype = numpy.float64  
index_dtype = numpy.uint16
bool_dtype = numpy.bool_

from westpa.binning.assign import UNKNOWN_INDEX as _UNKNOWN_INDEX
cdef index_t UNKNOWN_INDEX = _UNKNOWN_INDEX  

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef flux_assign(numpy.ndarray[weight_t, ndim=1] weights,
                  numpy.ndarray[index_t, ndim=1] init_assignments,
                  numpy.ndarray[index_t, ndim=1] final_assignments,
                  numpy.ndarray[weight_t, ndim=2] flux_matrix):
    cdef:
        Py_ssize_t m,n
        index_t i, j
    n = len(weights)
    for m from 0 <= m < n:
        i = init_assignments[m]
        j = final_assignments[m]
        flux_matrix[i,j] += weights[m]
    return

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef stats_process(numpy.ndarray[index_t, ndim=2] bin_assignments,
                    numpy.ndarray[weight_t, ndim=1] weights, 
                    numpy.ndarray[weight_t, ndim=1] fluxes, 
                    numpy.ndarray[weight_t, ndim=1] populations, 
                    numpy.ndarray[index_t, ndim=1] trans, 
                    numpy.ndarray[index_t, ndim=2] mask):
    cdef:
        Py_ssize_t i,k
        index_t ibin,fbin,nsegs,npts
    nsegs = bin_assignments.shape[0]
    npts = bin_assignments.shape[1]

    for i in xrange(0,npts - 1):
        for k in xrange(nsegs):
            ibin = bin_assignments[k,i]
            fbin = bin_assignments[k, i + 1]

            if mask[k, 0] == 1:
                continue

            w = weights[k]

            fluxes[ibin, fbin] += w
            trans[ibin, fbin] += 1
            populations[ibin] += w
    return
