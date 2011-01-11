__metaclass__ = type

from itertools import chain
import numpy

from wemd import Segment, Trajectory


class TrajTree:
    def __init__(self, data_manager, squeeze_data):
        self.data_manager = data_manager
        self.squeeze_data = squeeze_data
        self.roots = None
        
    def _seg_id(self, node):
        try:
            return node.seg_id
        except AttributeError:
            return long(node)
        
    def get_roots(self, max_iter = None):      
        if not max_iter:
            max_iter = self.data_manager.get_last_complete_iter()
              
        if self.roots is None:
            self.roots = self.data_manager.get_segments_by_parent_id(1, n_iter_upper = max_iter, p_parent_id = None)
                    
        return self.roots
    
    def get_children(self, node):
        seg_id = self._seg_id(node)
        max_iter = self.data_manager.get_last_complete_iter()          
        s = self.data_manager.get_segments_by_parent_id(1, n_iter_upper = max_iter, p_parent_id = seg_id)
        return s
    
    def _desc_traj(self, node, segs, max_iter, trajs):
        children = self.data_manager.get_segments_by_parent_id(1, n_iter_upper = max_iter, p_parent_id = node.seg_id)
        
        if not children or node.n_iter == max_iter:
            trajs.append(Trajectory(segs + [node], self.squeeze_data))
        else:
            for child in children:
                self._desc_traj(child, segs + [node], max_iter, trajs)
            
    
    def get_trajectories_recursive(self, max_iter):
        trajs = []
        for root in self.get_roots(max_iter):
            self._desc_traj(root, [], max_iter, trajs)
        return trajs
         
    def get_trajectories_iterative_reference(self, max_iter):
        trajs = []
        
        for root in self.get_roots(max_iter):
            state_stack = [{'node': root,
                           'children': self.data_manager.get_segments_by_parent_id(1, n_iter_upper = max_iter, p_parent_id = root.seg_id)
                           }]
            
                        
            while state_stack:
                state = state_stack.pop(-1)
                
                # We don't want to descend if we are already at max_iter
                # (end of the line)
                if state['node'].n_iter == max_iter:
                    #trajs.append(state['node'].seg_id)
                    trajs.append(Trajectory([si['node'] for si in state_stack] 
                                            + [state['node']],
                                            self.squeeze_data))
                else:                
                    while state['children']:
                        node = state['children'].pop(-1)                    
                        children = self.data_manager.get_segments_by_parent_id(1, n_iter_upper = max_iter, p_parent_id = node.seg_id)
                        if not children or node.n_iter == max_iter:
                            # end of the line
                            trajs.append(Trajectory([si['node'] for si in state_stack] 
                                                    + [state['node'], node],
                                                    self.squeeze_data))
                            #trajs.append(node.seg_id)
                        else:
                            state_stack.append(state)
                            state = {'node': node,
                                     'children': children}

        return trajs
    
    def trace_trajectories(self, max_iter, callable, get_ext_state = None,
                           set_ext_state = None):
        """
        Walk the trajectory tree depth-first, calling
          ``callable(segment, children, history)`` for each segment 
        visited. ``segment`` is the segment being visited, ``children`` is that
        segment's children, ``history`` is the chain of segments leading
        to ``segment`` (not including ``segment``).
        """
        
        if get_ext_state is None:
            if set_ext_state is not None:
                raise ValueError('external state getter required')
            else:
                get_ext_state = (lambda: None)
        if set_ext_state is None:
            if get_ext_state is not None:
                raise ValueError('external state setter required')
            else:
                set_ext_state = (lambda ext: None)
            
        
        # This will grow to contain the maximum trajectory length
        history_chunksize = 32
        history = numpy.empty((history_chunksize,), numpy.object_)
        
        #s = select([segments_table]).where((segments_table.c.p_parent_id == bindparam('p_parent_id')))
        
        for root in self.get_roots(max_iter):
            children = self.data_manager.get_segments_by_parent_id(1, n_iter_upper = max_iter, p_parent_id = root.seg_id)

            # Visit the root node of each tree unconditionally
            callable(root, children, [])

            
            state_stack = [{'node': root,
                            'children': children,
                            'len_history': 0,
                            'ext': get_ext_state()}]
            
            # Walk the tree, depth-first
            while state_stack:
                state = state_stack.pop(-1)
                
                node = state['node']
                children = state['children']
                len_history = state['len_history']
                set_ext_state(state['ext'])
                
                # Descend as far as we can
                while node.n_iter < max_iter and len(children):
                    # Save current state before descending
                    state_stack.append({'node': node,
                                        'children': children,
                                        'len_history': len_history,
                                        'ext': get_ext_state()})
                    
                    # Add an item to the historical record
                    if len_history >= history.shape[0]: 
                        history.resize((history.shape[0] + history_chunksize,))
                    history[len_history] = node
                    len_history += 1
                    
                    node = children.pop(-1)
                    children = self.data_manager.get_segments_by_parent_id(1, n_iter_upper = max_iter, p_parent_id = node.seg_id)
                    
                    # Visit the new node as we descend
                    callable(node, children, history[:len_history])
                    
    def get_trajectories(self, max_iter):
        trajs = []
        def append_traj(node, children, history):
            if node.n_iter == max_iter or not children:
                trajs.append(Trajectory(chain(history, [node]), self.squeeze_data))
        self.trace_trajectories(max_iter, append_traj)
        return trajs

        