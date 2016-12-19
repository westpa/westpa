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

from __future__ import division, print_function; __metaclass__ = type
from itertools import izip
from collections import deque

cimport cython
import numpy
cimport numpy
from numpy cimport uint32_t, uint64_t, int64_t, float64_t


_tt_dtype = numpy.dtype([('n_iter', numpy.uint32),
                         ('seg_id', numpy.int64),
                         ('parent_id', numpy.int64),
                         ('parent_offset', numpy.int64),  # offset of parent segment into this table
                         ('weight', numpy.float64),       # weight of this segment
                         ])

cdef packed struct _tt_rec:
    uint32_t  n_iter
    int64_t   seg_id
    int64_t   parent_id
    int64_t   parent_offset
    float64_t weight

cdef int64_t NO_PARENT = -1  #indicates that no parent is contained in the table; segment is a root
cdef uint32_t CT_CHUNKSIZE = 512
cdef uint32_t NODE_CHUNKSIZE = 512

cdef class _trajtree_base:
    cdef public object trajtable    # numpy.ndarray[_tt_rec, ndim=1]
    cdef public object childtable   # numpy.ndarray[object, ndim=1]
    cdef public object iter_offsets # dict
    
    #@cython.boundscheck(False)
    #@cython.wraparound(False)
    cpdef _build_table(self, segsel, data_manager):
    
        cdef uint32_t n_iter
        cdef uint32_t start_iter = segsel.start_iter
        cdef uint32_t stop_iter = segsel.stop_iter
        cdef uint64_t seg_count = len(segsel)
        
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = numpy.empty((seg_count,), dtype=_tt_dtype)
        cdef numpy.ndarray[object, ndim=1]  childtable = numpy.empty((seg_count,), dtype=numpy.object_)
        cdef numpy.ndarray[int64_t, ndim=1] seg_ids = None
        cdef numpy.ndarray[int64_t, ndim=1] all_parent_ids = None
        cdef numpy.ndarray[float64_t, ndim=1] weights = None 
        
        cdef uint64_t tt_offset = 0
        cdef int64_t seg_id = 0
        cdef int64_t parent_id = 0
        cdef int64_t parent_offset = NO_PARENT
        
        last_iter_indices = {} # mapping of seg_id -> table row for that seg_id in previous iteration
        iter_offsets = {}
        for n_iter in range(start_iter, stop_iter):
            iter_offsets[n_iter] = tt_offset
            this_iter_indices = {}
            
            seg_ids = numpy.fromiter(segsel.from_iter(n_iter), dtype=numpy.int64)
            iter_group = data_manager.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            weights = seg_index['weight']
            all_parent_ids = numpy.array(data_manager.get_all_parent_ids(n_iter), numpy.int64)
            
            for seg_id in seg_ids:
                this_iter_indices[seg_id] = tt_offset
                
                # Record information about this segment
                trajtable[tt_offset].n_iter = n_iter
                trajtable[tt_offset].seg_id = seg_id
                trajtable[tt_offset].weight = weights[seg_id]
                trajtable[tt_offset].parent_id = all_parent_ids[seg_id]
                childtable[tt_offset] = None

                # Now record information about this segment's parent
                if n_iter == start_iter or trajtable[tt_offset].parent_id < 0:
                    trajtable[tt_offset].parent_offset = NO_PARENT
                else:
                    # where in the table is the parent?
                    # note that this will raise KeyError if a segment required for proper connectivity is
                    # missing from the input segment set
                    parent_offset = last_iter_indices[trajtable[tt_offset].parent_id]
                    #assert parent_offset != -1 and parent_offset < tt_offset
                    #assert trajtable[parent_offset].seg_id == trajtable[tt_offset].parent_id
                    #assert trajtable[parent_offset].n_iter == n_iter-1
                    trajtable[tt_offset].parent_offset = parent_offset
                    
                    # all first-iteration segments are roots for the purposes of analysis
                    # as are true roots (from recycling) 
                    if childtable[parent_offset] is None:
                        childtable[parent_offset] = [tt_offset]
                    else:
                        childtable[parent_offset].append(tt_offset)
                tt_offset += 1
                
            del last_iter_indices, seg_index                                 
            last_iter_indices = this_iter_indices
        
        for tt_offset in range(len(trajtable)):
            if childtable[tt_offset] is not None:
                childtable[tt_offset] = numpy.array(childtable[tt_offset], dtype=numpy.int64)

        self.childtable = childtable
        self.trajtable = trajtable    
        self.iter_offsets = iter_offsets
        
    cpdef uint64_t count_roots(self):
        cdef uint64_t i=0, root_count=0
        
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = self.trajtable
        
        for i in range(0, len(trajtable)):
            if trajtable[i].parent_offset == NO_PARENT:
                root_count += 1
        return root_count
    
    cpdef uint64_t count_leaves(self):
        cdef uint64_t i=0, leaf_count=0
        
        cdef numpy.ndarray[object, ndim=1] childtable = self.childtable
        for i in range(0, len(childtable)):
            if childtable[i] is None:
                leaf_count += 1
        return leaf_count

    cpdef get_child_indices(self, int64_t index):
        child_indices = self.childtable[index]
        if child_indices is None:
            return []
        else:
            return child_indices
