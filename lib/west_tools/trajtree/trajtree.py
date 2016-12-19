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

from collections import namedtuple

import numpy
trajnode = namedtuple('trajnode', ('n_iter', 'seg_id'))
from collections import deque
import westpa
from westtools.tool_classes.selected_segs import AllSegmentSelection
import _trajtree
from _trajtree import _trajtree_base #@UnresolvedImport

class TrajTreeSet(_trajtree_base):    
    def __init__(self, segsel = None, data_manager = None):
        self.data_manager = data_manager or westpa.rc.get_data_manager()
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

                state_stack.append({'subtrees': subtrees,
                                    'vstate': get_visitor_state() if get_visitor_state else None})
                
                subtrees = deque(self.get_child_indices(index))
                
                n_visits += 1
                try:
                    visit(node['n_iter'], node['seg_id'], node['weight'], has_children = (len(subtrees) > 0), *vargs, **vkwargs)
                except StopIteration:
                    subtrees = deque()
                    continue # to next sibling
        
        return n_visits


class FakeTrajTreeSet(TrajTreeSet):
    def __init__(self):

#_tt_dtype = numpy.dtype([('n_iter', numpy.uint32),
#                         ('seg_id', numpy.int64),
#                         ('parent_id', numpy.int64),
#                         ('parent_offset', numpy.int64),  # offset of parent segment into this table
#                         ('weight', numpy.float64)])       # weight of this segment

        
        self.trajtable = numpy.array([(1, 1,-1,-1,1.0), #0
                                      (1,11,-1,-1,1.0), #1
                                      (2, 2, 1, 0,1.0), #2
                                      (2, 3, 1, 0,1.0), #3
                                      (2, 4, 1, 0,1.0), #4
                                      (2,12,11, 1,1.0), #5
                                      (2,13,11, 1,1.0), #6
                                      (3, 5, 2, 2,1.0), #7
                                      (3, 6, 3, 3,1.0), #8
                                      (3, 7, 4, 4,1.0), #9
                                      (3, 8, 4, 4,1.0), #10
                                      (3,14,12, 5,1.0), #11
                                      (3,15,12, 5,1.0), #12
                                      (4, 9, 5, 7,1.0), #13
                                      (4,10, 5, 7,1.0), #14
                                      ], dtype=_trajtree._tt_dtype)

        empty_array = numpy.array([])
        self.childtable = numpy.array([numpy.array([2, 3, 4]), # segment 1
                                       numpy.array([5, 6]),    # segment 11
                                       numpy.array([7]),       # segment 2
                                       numpy.array([8]),       # segment 3
                                       numpy.array([9, 10]),   # segment 4
                                       numpy.array([11, 12]),  # segment 12
                                       empty_array,            # segment 13
                                       numpy.array([13, 14]),  # segment 5
                                       empty_array,            # segment 6
                                       empty_array,            # segment 7
                                       empty_array,            # segment 8
                                       empty_array,            # segment 14
                                       empty_array,            # segment 15
                                       empty_array,            # segment 9
                                       empty_array,            # segment 10
                                       ], dtype=numpy.object_)
        self.iter_offsets = {1: 0,
                             2: 2,
                             3: 7,
                             4: 13}
        
