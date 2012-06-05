from collections import namedtuple

import numpy
trajnode = namedtuple('trajnode', ('n_iter', 'seg_id'))
from collections import deque
import wemd
from wt2.tool_classes.selected_segs import AllSegmentSelection
import _trajtree
from _trajtree import _trajtree_base
from _trajtree import node_dtype

class TrajTreeSet(_trajtree_base):    
    def __init__(self, segsel = None, data_manager = None):
        self.data_manager = data_manager or wemd.rc.get_data_manager()
        self.segsel = segsel or AllSegmentSelection(data_manager = self.data_manager)
        self._build_table()
        
    def _build_table(self):
        _trajtree_base._build_table(self, self.segsel, self.data_manager)    
                
    def __len__(self):
        return len(self.trajtable)
    
    def count_roots(self):
        return self._count_roots()
    
    def count_leaves(self):
        return self._count_leaves()
    
    def get_roots(self):
        return [trajnode(root['n_iter'], root['seg_id']) for root in self._get_roots()]
    
    def get_child_ids(self, n_iter, seg_id):
        tt_offset = self.iter_offsets[n_iter]
        trajtable = self.trajtable
        
        while trajtable[tt_offset]['seg_id'] != seg_id:
            tt_offset += 1
            
        offset = trajtable[tt_offset]['children_offset']
        n_children = trajtable[tt_offset]['n_children']
        if n_children == 0:
            return []
        
        return trajtable['seg_id'][self.childtable[offset:offset+n_children]]
        
    
#    def count_roots(self):
#        return (self.trajtable['parent_offset'] == -1).sum()
#    
#    def count_leaves(self):
#        return (self.trajtable['n_children'] == 0).sum()

    def trace_trajectories(self, visit, get_visitor_state = None, set_visitor_state = None, vargs=None, vkwargs=None):
        
        if (get_visitor_state or set_visitor_state) and not (get_visitor_state and set_visitor_state):
            raise ValueError('either both or neither of get_visitor_state and set_visitor_state must be specified')
        
        vargs = vargs or ()
        vkwargs = vkwargs or {}
        n_visits = 0
        
        roots = deque(self.get_roots())
        print('Examining {:d} roots'.format(len(roots)))
        state_stack = deque([{'subtrees': roots,
                              'vstate': get_visitor_state() if get_visitor_state else None}])
        while state_stack:
            state = state_stack.pop()
            subtrees = state['subtrees']
            if set_visitor_state:
                set_visitor_state(state['vstate'])
            
            while subtrees:
                node = subtrees.popleft()
                n_visits += 1
                visit(node.n_iter, node.seg_id, *vargs, **vkwargs)
                
                state_stack.append({'subtrees': subtrees,
                                    'vstate': get_visitor_state() if get_visitor_state else None})
                
                subtrees = deque(trajnode(node.n_iter+1, child_id) for child_id in self.get_child_ids(node.n_iter, node.seg_id))
        
        return n_visits
            
            # next state level

#_tt_dtype = numpy.dtype([('n_iter', numpy.uint32),
#                         ('seg_id', numpy.int64),
#                         ('parent_id', numpy.int64),
#                         ('parent_offset', numpy.int64),  # offset of parent segment into this table
#                         ('n_children', numpy.uint32),     # number of children
#                         ('children_offset', numpy.int64),# offset of child segments into child table
#                         ('weight', numpy.float64)])       # weight of this segment


class FakeTrajTreeSet(TrajTreeSet):
    def __init__(self):
        
        self.trajtable = numpy.array([(1, 1,-1,-1,3, 0,1.0), #0
                                      (1,11,-1,-1,2, 3,1.0), #1
                                      (2, 2, 1, 0,1, 5,1.0), #2
                                      (2, 3, 1, 0,1, 6,1.0), #3
                                      (2, 4, 1, 0,2, 7,1.0), #4
                                      (2,12,11, 1,2, 9,1.0), #5
                                      (2,13,11, 1,0, 0,1.0), #6
                                      (3, 5, 2, 2,2,11,1.0), #7
                                      (3, 6, 3, 3,0, 0,1.0), #8
                                      (3, 7, 4, 4,0, 0,1.0), #9
                                      (3, 8, 4, 4,0, 0,1.0), #10
                                      (3,14,12, 5,0, 0,1.0), #11
                                      (3,15,12, 5,0, 0,1.0), #12
                                      (4, 9, 5, 7,0, 0,1.0), #13
                                      (4,10, 5, 7,0, 0,1.0), #14
                                      ], dtype=_trajtree._tt_dtype)

        self.childtable = numpy.array([2, 3, 4, # 0: segment 1
                                       5, 6,    # 3: segment 11
                                       7,       # 5: segment 2
                                       8,       # 6: segment 3
                                       9, 10,   # 7: segment 4
                                       11, 12,  # 9: segment 12
                                       13, 14,  #11: segment 5
                                       ], dtype=numpy.uint64)
        self.iter_offsets = {1: 0,
                             2: 2,
                             3: 7,
                             4: 13}
        
