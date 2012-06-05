from __future__ import division, print_function; __metaclass__ = type
from itertools import izip

cimport cython
import numpy
cimport numpy
from numpy cimport uint32_t, uint64_t, int64_t, float64_t


_tt_dtype = numpy.dtype([('n_iter', numpy.uint32),
                         ('seg_id', numpy.int64),
                         ('parent_id', numpy.int64),
                         ('parent_offset', numpy.int64),  # offset of parent segment into this table
                         ('n_children', numpy.uint32),     # number of children
                         ('children_offset', numpy.int64),# offset of child segments into child table
                         ('weight', numpy.float64)])       # weight of this segment

cdef packed struct _tt_rec:
    uint32_t  n_iter
    int64_t   seg_id
    int64_t   parent_id
    int64_t   parent_offset
    uint32_t  n_children
    int64_t   children_offset
    float64_t weight
    
node_dtype = numpy.dtype([('n_iter', numpy.uint32),
                          ('seg_id', numpy.int64)])

cdef packed struct _node_rec:
    uint32_t n_iter
    int64_t  seg_id     

cdef int64_t NO_PARENT = -1  #indicates that no parent is contained in the table; segment is a root
cdef uint32_t CT_CHUNKSIZE = 512
cdef uint32_t NODE_CHUNKSIZE = 512

cdef class _trajtree_base:
    cdef public object trajtable    # numpy.ndarray[_tt_rec, ndim=1]
    cdef public object childtable   # numpy.ndarray[numpy.uint64_t, ndim=1]
    cdef public object iter_offsets # dict
    
    #@cython.boundscheck(False)
    #@cython.wraparound(False)
    cpdef _build_table(self, segsel, data_manager):
    
        cdef uint32_t n_iter
        cdef uint32_t start_iter = segsel.start_iter
        cdef uint32_t stop_iter = segsel.stop_iter
        cdef uint64_t seg_count = len(segsel)
        
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = numpy.empty((seg_count,), dtype=_tt_dtype)
        childtable = numpy.empty((CT_CHUNKSIZE,), dtype=numpy.uint64)
        cdef numpy.ndarray[int64_t, ndim=1] seg_ids = None
        #cdef numpy.ndarray[int64_t, ndim=1] all_parent_ids = None
        cdef numpy.ndarray[float64_t, ndim=1] weights = None 
        
        cdef uint64_t len_childtable = 0
        cdef uint64_t tt_offset = 0
        cdef int64_t seg_id = 0
        cdef int64_t parent_id = 0
        cdef int64_t parent_offset = NO_PARENT
        
        last_iter_indices = {} # mapping of seg_id -> table row for that seg_id in previous iteration
        iter_offsets = {}
        for n_iter in range(start_iter, stop_iter):
            iter_offsets[n_iter] = tt_offset
            this_iter_indices = {}
            
            seg_ids = numpy.array(segsel.from_iter(n_iter), dtype=numpy.int64)
            iter_group = data_manager.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            weights = seg_index['weight']
            #all_parent_ids = numpy.array(data_manager.get_all_parent_ids(n_iter), dtype=numpy.int64)
            all_parent_ids = data_manager.get_all_parent_ids(n_iter)
            
            for seg_id in seg_ids:
                this_iter_indices[seg_id] = tt_offset
                
                # Record information about this segment
                trajtable[tt_offset].n_iter = n_iter
                trajtable[tt_offset].seg_id = seg_id
                trajtable[tt_offset].weight = weights[seg_id]
                trajtable[tt_offset].n_children = 0
                trajtable[tt_offset].children_offset = 0
                
                # Now record information about this segment's parent
                if n_iter == start_iter:
                    # all first-iteration segments are roots for the purposes of analysis
                    trajtable[tt_offset].parent_offset = NO_PARENT
                else:
                    # who is our parent?
                    parent_id = all_parent_ids[seg_id]
                    trajtable[tt_offset].parent_id = parent_id
                    if parent_id < 0:
                        trajtable[tt_offset].parent_offset = NO_PARENT
                    else:
                        trajtable[tt_offset].parent_offset = last_iter_indices[parent_id]
                    
                # And now incorporate information about whose child this segment is
                parent_offset = trajtable[tt_offset].parent_offset
                if parent_offset != NO_PARENT:
                    # need to add a child to the list
                    if len_childtable == len(childtable):
                        childtable.resize((len_childtable+CT_CHUNKSIZE,))
                    childtable[len_childtable] = tt_offset
                    
                    if trajtable[parent_offset].n_children == 0:
                        trajtable[parent_offset].n_children = 1
                        trajtable[parent_offset].children_offset = len_childtable
                    else:
                        trajtable[parent_offset].n_children += 1 
                    
                    len_childtable+=1
                
                tt_offset += 1

            del last_iter_indices                                 
            last_iter_indices = this_iter_indices
            
        childtable.resize((len(childtable),))
        self.childtable = childtable
        self.trajtable = trajtable    
        self.iter_offsets = iter_offsets
        
    cpdef uint64_t _count_roots(self):
        cdef uint64_t i=0, root_count=0
        
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = self.trajtable
        
        for i in range(0, len(trajtable)):
            if trajtable[i].parent_offset == NO_PARENT:
                root_count += 1
                
        return root_count
    
    cpdef uint64_t _count_leaves(self):
        cdef uint64_t i=0, leaf_count=0
        
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = self.trajtable
        
        for i in range(0, len(trajtable)):
            if trajtable[i].n_children == 0:
                leaf_count += 1
                
        return leaf_count
    
    cpdef object _get_roots(self):
        cdef int64_t i=0, n_roots=0
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = self.trajtable        
        roots = numpy.empty((NODE_CHUNKSIZE,), dtype=node_dtype)
        
        for i in range(0, len(trajtable)):
            if trajtable[i].parent_offset == NO_PARENT:
                if n_roots == len(roots):
                    roots.resize((n_roots+NODE_CHUNKSIZE,))
                roots[n_roots] = trajtable[i].n_iter, trajtable[i].seg_id
                n_roots += 1
        
        roots.resize((n_roots,))
        #roots = numpy.resize(roots , (n_roots,))
        return roots
            
    cpdef object _get_child_ids(self, n_iter, seg_id):
    
        cdef uint32_t _n_iter = n_iter
        cdef int64_t  _seg_id = seg_id
        cdef uint64_t tt_offset, offset, length
        cdef numpy.ndarray[_tt_rec, ndim=1] trajtable = self.trajtable
        
        tt_offset = self.iter_offsets[_n_iter]
        
        while trajtable[tt_offset].seg_id != _seg_id:
            tt_offset += 1
            assert tt_offset <= len(trajtable)
            
        offset = trajtable[tt_offset].children_offset
        length = trajtable[tt_offset].n_children

        sl = trajtable.take(self.childtable[offset:offset+length], axis=-1)
        print(n_iter, seg_id, sl)
        for row in sl:
            assert row['parent_id'] == seg_id
            assert row['n_iter'] == n_iter+1
        return sl           
        #return self.childtable[offset:offset+length]
        
                
                
                
            
            
        
        