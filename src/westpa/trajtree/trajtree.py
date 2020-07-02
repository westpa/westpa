import collections

import numpy as np

import westpa
from westpa.tools.selected_segs import AllSegmentSelection
from . import _trajtree
from ._trajtree import _trajtree_base  # @UnresolvedImport


trajnode = collections.namedtuple('trajnode', ('n_iter', 'seg_id'))


class TrajTreeSet(_trajtree_base):
    def __init__(self, segsel=None, data_manager=None):
        self.data_manager = data_manager or westpa.rc.get_data_manager()
        self.segsel = segsel or AllSegmentSelection(data_manager=self.data_manager)
        self._build_table(self.segsel, self.data_manager)

    def __len__(self):
        return len(self.trajtable)

    def get_roots(self):
        return self.trajtable[self.trajtable['parent_offset'] == -1]
        # return [trajnode(root['n_iter'], root['seg_id']) for root in self._get_roots()]

    def get_root_indices(self):
        return np.squeeze(np.argwhere(self.trajtable['parent_offset'] == -1))

    def trace_trajectories(self, visit, get_visitor_state=None, set_visitor_state=None, vargs=None, vkwargs=None):

        if (get_visitor_state or set_visitor_state) and not (get_visitor_state and set_visitor_state):
            raise ValueError('either both or neither of get_visitor_state and set_visitor_state must be specified')

        vargs = vargs or ()
        vkwargs = vkwargs or {}
        n_visits = 0

        trajtable = self.trajtable

        roots = collections.deque(self.get_root_indices())
        print('Examining {:d} roots'.format(len(roots)))
        state_stack = collections.deque([{'subtrees': roots, 'vstate': get_visitor_state() if get_visitor_state else None}])
        while state_stack:
            state = state_stack.pop()
            subtrees = state['subtrees']
            if set_visitor_state:
                set_visitor_state(state['vstate'])

            while subtrees:
                index = subtrees.popleft()
                node = trajtable[index]

                state_stack.append({'subtrees': subtrees, 'vstate': get_visitor_state() if get_visitor_state else None})

                subtrees = collections.deque(self.get_child_indices(index))

                n_visits += 1
                try:
                    visit(node['n_iter'], node['seg_id'], node['weight'], has_children=(len(subtrees) > 0), *vargs, **vkwargs)
                except StopIteration:
                    subtrees = collections.deque()
                    continue  # to next sibling

        return n_visits


class FakeTrajTreeSet(TrajTreeSet):
    def __init__(self):

        # _tt_dtype = np.dtype([('n_iter', np.uint32),
        #                         ('seg_id', np.int64),
        #                         ('parent_id', np.int64),
        #                         ('parent_offset', np.int64),  # offset of parent segment into this table
        #                         ('weight', np.float64)])       # weight of this segment

        self.trajtable = np.array(
            [
                (1, 1, -1, -1, 1.0),  # 0
                (1, 11, -1, -1, 1.0),  # 1
                (2, 2, 1, 0, 1.0),  # 2
                (2, 3, 1, 0, 1.0),  # 3
                (2, 4, 1, 0, 1.0),  # 4
                (2, 12, 11, 1, 1.0),  # 5
                (2, 13, 11, 1, 1.0),  # 6
                (3, 5, 2, 2, 1.0),  # 7
                (3, 6, 3, 3, 1.0),  # 8
                (3, 7, 4, 4, 1.0),  # 9
                (3, 8, 4, 4, 1.0),  # 10
                (3, 14, 12, 5, 1.0),  # 11
                (3, 15, 12, 5, 1.0),  # 12
                (4, 9, 5, 7, 1.0),  # 13
                (4, 10, 5, 7, 1.0),  # 14
            ],
            dtype=_trajtree._tt_dtype,
        )

        empty_array = np.array([])
        self.childtable = np.array(
            [
                np.array([2, 3, 4]),  # segment 1
                np.array([5, 6]),  # segment 11
                np.array([7]),  # segment 2
                np.array([8]),  # segment 3
                np.array([9, 10]),  # segment 4
                np.array([11, 12]),  # segment 12
                empty_array,  # segment 13
                np.array([13, 14]),  # segment 5
                empty_array,  # segment 6
                empty_array,  # segment 7
                empty_array,  # segment 8
                empty_array,  # segment 14
                empty_array,  # segment 15
                empty_array,  # segment 9
                empty_array,  # segment 10
            ],
            dtype=np.object_,
        )
        self.iter_offsets = {1: 0, 2: 2, 3: 7, 4: 13}
