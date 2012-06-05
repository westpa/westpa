from __future__ import division, print_function; __metaclass__ = type
from itertools import izip

cimport cython
import numpy
cimport numpy
from numpy cimport uint32_t, uint64_t, int64_t, float64_t


_tt_dtype = numpy.dtype([('n_iter', numpy.uint32),
                         ('seg_id', numpy.int64),
                         ('parent_offset', numpy.int64),  # offset of parent segment into this table
                         ('n_children', numpy.uint32),     # number of children
                         ('children_offset', numpy.int64),# offset of child segments into child table
                         ('weight', numpy.float64)])       # weight of this segment

cdef packed struct _tt_rec:
    uint32_t  n_iter
    int64_t  seg_id
    int64_t  parent_offset
    uint32_t  n_children
    int64_t  children_offset
    float64_t weight

cdef int64_t NO_PARENT = -1  #indicates that no parent is contained in the table; segment is a root
cdef uint32_t CT_CHUNKSIZE = 512

cdef class _trajtree_base:
    cdef readonly object trajtable    # numpy.ndarray[_tt_rec, ndim=1]
    cdef readonly object childtable   # numpy.ndarray[numpy.uint64_t, ndim=1]
    cdef readonly object iter_offsets # dict
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _build_table(self, segsel, data_manager):
    
        cdef uint32_t n_iter
        cdef uint32_t start_iter = segsel.start_iter
        cdef uint64_t seg_count = len(segsel)
        
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = numpy.empty((seg_count,), dtype=_tt_dtype)
        childtable = numpy.empty((CT_CHUNKSIZE,), dtype=numpy.uint64)
        cdef numpy.ndarray[int64_t, ndim=1] all_parent_ids = None
        cdef numpy.ndarray[float64_t, ndim=1] weights = None 
        
        cdef uint64_t len_childtable = 0
        cdef uint64_t tt_offset = 0
        cdef int64_t seg_id = 0
        cdef int64_t parent_id = 0
        cdef int64_t parent_row = NO_PARENT
        
        last_iter_indices = {} # mapping of seg_id -> table row for that seg_id in previous iteration
        self.iter_offsets = {}
        for n_iter in xrange(start_iter, segsel.stop_iter):
            self.iter_offsets[n_iter] = tt_offset
            this_iter_indices = {}
            
            seg_ids = segsel.from_iter(n_iter)
            iter_group = data_manager.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            weights = seg_index['weight']
            all_parent_ids = numpy.array(data_manager.get_all_parent_ids(n_iter), dtype=numpy.int64)
            
            for seg_id in seg_ids:
                this_iter_indices[seg_id] = tt_offset
                trajtable[tt_offset].n_iter = n_iter
                trajtable[tt_offset].seg_id = seg_id
                trajtable[tt_offset].weight = weights[seg_id]
                trajtable[tt_offset].n_children = 0
                trajtable[tt_offset].children_offset = 0
            
                if n_iter == start_iter:
                    parent_row = NO_PARENT
                else:
                    parent_id = all_parent_ids[seg_id]
                    
                    if parent_id < 0:
                        parent_row = NO_PARENT
                    else:                        
                        # this will raise KeyError if segsel is missing segments along a trajectory
                        parent_row = last_iter_indices[parent_id]
                        assert parent_row != NO_PARENT
        
                trajtable[tt_offset].parent_offset = parent_row
                
                if parent_row != NO_PARENT:
                    if trajtable[parent_row].n_children == 0:
                        # First child
                        trajtable[parent_row].n_children = 1
                        trajtable[parent_row].children_offset = len_childtable
                    else:
                        trajtable[parent_row].n_children += 1
                    
                    if len_childtable == len(childtable):
                        childtable.resize((len(childtable)+CT_CHUNKSIZE,))
                        childtable[len_childtable] = tt_offset
                    len_childtable += 1
                
                tt_offset += 1
                                 
            last_iter_indices = this_iter_indices
            
        childtable.resize((len(childtable),))
        self.childtable = childtable
        self.trajtable = trajtable    
    
            
            
        
                
                
                
            
            
        
        