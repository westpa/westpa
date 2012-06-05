from collections import namedtuple

import numpy
trajnode = namedtuple('trajnode', ('n_iter', 'seg_id'))
from collections import deque
import wemd
from wt2.tool_classes.selected_segs import AllSegmentSelection
import _trajtree
from _trajtree import _trajtree_base

class TrajTreeSet(_trajtree_base):    
    def __init__(self, segsel = None, data_manager = None):
        self.data_manager = data_manager or wemd.rc.get_data_manager()
        self.segsel = segsel or AllSegmentSelection(data_manager = self.data_manager)
        self._build_table(self.segsel, self.data_manager)
                        
    def __len__(self):
        return len(self.trajtable)
        
    def get_roots(self):
        return self.trajtable[self.trajtable['parent_offset'] == -1]
        #return [trajnode(root['n_iter'], root['seg_id']) for root in self._get_roots()]
        
    def get_root_indices(self):
        return numpy.squeeze(numpy.argwhere(self.trajtable['parent_offset'] == -1))
            
    def trace_trajectories(self, visit, get_visitor_state = None, set_visitor_state = None, vargs=None, vkwargs=None):
        
        if (get_visitor_state or set_visitor_state) and not (get_visitor_state and set_visitor_state):
            raise ValueError('either both or neither of get_visitor_state and set_visitor_state must be specified')
        
        vargs = vargs or ()
        vkwargs = vkwargs or {}
        n_visits = 0
        
        trajtable = self.trajtable
        
        roots = deque(self.get_root_indices())
        print('Examining {:d} roots'.format(len(roots)))
        state_stack = deque([{'subtrees': roots,
                              'vstate': get_visitor_state() if get_visitor_state else None}])
        while state_stack:
            state = state_stack.pop()
            subtrees = state['subtrees']
            if set_visitor_state:
                set_visitor_state(state['vstate'])
            
            while subtrees:
                index = subtrees.popleft()
                node = trajtable[index]
                n_visits += 1
                visit(node['n_iter'], node['seg_id'], *vargs, **vkwargs)
                
                state_stack.append({'subtrees': subtrees,
                                    'vstate': get_visitor_state() if get_visitor_state else None})
                
                #subtrees = deque(trajnode(node.n_iter+1, child_id) for child_id in self.get_child_ids(node.n_iter, node.seg_id))
                subtrees = deque(self.get_child_indices(index))
        
        return n_visits
            



class FakeTrajTreeSet(TrajTreeSet):
    def __init__(self):

#_tt_dtype = numpy.dtype([('n_iter', numpy.uint32),
#                         ('seg_id', numpy.int64),
#                         ('parent_id', numpy.int64),
#                         ('parent_offset', numpy.int64),  # offset of parent segment into this table
#                         ('n_children', numpy.uint32),     # number of children
#                         ('children_offset', numpy.int64),# offset of child segments into child table
#                         ('weight', numpy.float64)])       # weight of this segment

        
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
        
