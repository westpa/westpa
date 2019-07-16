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


import logging

log = logging.getLogger(__name__)

import numpy

class TrajWalker:
    """A class to perform analysis by walking the trajectory tree.  A stack is used rather than recursion, or else
    the highest number of iterations capable of being considered would be the same as the Python recursion limit.
    """
    
    def __init__(self, data_reader, history_chunksize = 100):
        self.data_reader = data_reader        
        self.history_chunksize = history_chunksize
        self.n_segs_visited = 0        

    # TrajTree.count_segs_in_range() is now DataReader.total_segs_in_range()                

    def trace_to_root(self, n_iter, seg_id):
        '''Trace the given segment back to its starting point, returning a list of Segment
        objects describing the entire trajectory.'''
        
        segments = []
        segment = self.data_reader.get_segments_by_id(n_iter, [seg_id])[0]
        segments.append(segment)
        while segment.p_parent_id >= 0:
            segment = self.data_reader.get_segments_by_id(segment.n_iter-1, [segment.p_parent_id])[0]
            segments.append(segment)
        return list(reversed(segments))
                
    def get_trajectory_roots(self, first_iter, last_iter, include_pcoords = True):
        '''Get segments which start new trajectories.  If min_iter or max_iter is specified, restrict the
        set of iterations within which the search is conducted.'''
                        
        roots = []
        for n_iter in range(first_iter, last_iter+1):
            seg_ids = self.data_reader.get_created_seg_ids(n_iter)
            segments = self.data_reader.get_segments_by_id(n_iter, seg_ids, include_pcoords = include_pcoords)
            roots.extend(segments)
        return roots
    
    def get_initial_nodes(self, first_iter, last_iter, include_pcoords = True):
        '''Get segments with which to begin a tree walk -- those alive or created within [first_iter,last_iter].'''
                    
        root_ids = dict()
        
        # All trajectories alive or newly created in first_iter are initial nodes
        root_ids[first_iter] = set(self.data_reader.get_seg_ids(first_iter))
        
        # Find trajectories created in [first_iter, last_iter]
        for n_iter in range(first_iter, last_iter+1):
            seg_ids = self.data_reader.get_created_seg_ids(n_iter)
            try:
                root_ids[n_iter].update(seg_ids)
            except KeyError:
                root_ids[n_iter] = set(seg_ids)
                
        # Convert to Segment objects
        segments = []
        for (n_iter, id_set) in root_ids.items():
            segments.extend(self.data_reader.get_segments_by_id(n_iter, id_set, include_pcoords = include_pcoords))
        return segments
    
    def trace_trajectories(self, first_iter, last_iter, callable, include_pcoords = True, cargs = None, ckwargs = None,  
                           get_state = None, set_state = None):
        """
        Walk the trajectory tree depth-first, calling
          ``callable(segment, children, history, *cargs, **ckwargs)`` for each segment 
        visited. ``segment`` is the segment being visited, ``children`` is that
        segment's children, ``history`` is the chain of segments leading
        to ``segment`` (not including ``segment``). get_state and set_state are
        used to record and reset, respectively, any state specific to 
        ``callable`` when a new branch is traversed.
        """
        
        cargs = cargs or tuple()
        ckwargs = ckwargs or dict()
        
        # Either both or neither of external state getter/setter required
        if (get_state or set_state) and not (get_state and set_state):
            raise ValueError('either both or neither of get_state/set_state must be specified')
            
        # This will grow to contain the maximum trajectory length
        history = numpy.empty((self.history_chunksize,), numpy.object_)
        roots = self.get_initial_nodes(first_iter, last_iter, include_pcoords)

        for root in roots:
            children = self.data_reader.get_children(root, include_pcoords)

            # Visit the root node of each tree unconditionally
            callable(root, children, [], *cargs, **ckwargs)
            self.n_segs_visited += 1

            state_stack = [{'node': root,
                            'children': children,
                            'len_history': 0,
                            'ext': get_state() if get_state else None}]
            
            # Walk the tree, depth-first
            while state_stack:
                state = state_stack.pop(-1)
                
                node = state['node']
                children = state['children']
                len_history = state['len_history']
                if set_state: set_state(state['ext'])
                
                # Descend as far as we can
                while node.n_iter < last_iter and len(children):
                    # Save current state before descending
                    state_stack.append({'node': node,
                                        'children': children,
                                        'len_history': len_history,
                                        'ext': get_state() if get_state else None})
                    
                    # Add an item to the historical record
                    if len_history >= history.shape[0]: 
                        history.resize((history.shape[0] + self.history_chunksize,))
                    history[len_history] = node
                    len_history += 1
                    
                    node = children.pop(-1)
                    children = self.data_reader.get_children(node, include_pcoords)

                    # Visit the new node as we descend
                    callable(node, children, history[:len_history], *cargs, **ckwargs)
                    self.n_segs_visited += 1
                    