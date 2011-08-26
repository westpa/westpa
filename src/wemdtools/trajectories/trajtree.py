from __future__ import print_function, division; __metaclass__=type
from itertools import imap
import numpy
from wemd import Segment
from wemdtools.data_manager import CachingWEMDDataReader

import logging
log = logging.getLogger(__name__)

class TrajTree:
    def __init__(self, data_manager, include_pcoords = False, cache_pcoords = False):
        if not isinstance(data_manager,CachingWEMDDataReader):
            self.data_manager = CachingWEMDDataReader(data_manager, cache_pcoords)
        else:
            self.data_manager = data_manager
        self.data_manager.include_pcoords = include_pcoords
        
        self.history_chunksize = 32        
        self.segments_to_visit = 0
        self.segments_visited = 0
        
    def count_segs_in_range(self, min_iter = None, max_iter = None):
        min_iter = min_iter or 1
        if max_iter is None or max_iter >= self.data_manager.current_iteration:
            max_iter = self.data_manager.current_iteration - 1 
        
        self.segments_to_visit = 0
        for n_iter in xrange(min_iter, max_iter+1):
            seg_index = self.data_manager.get_seg_index(n_iter)
            self.segments_to_visit += len(seg_index)        
                
    def get_roots(self, min_iter = None, max_iter = None):
        '''Get segments which start new trajectories.  If min_iter or max_iter is specified, restrict the
        set of iterations within which the search is conducted.'''
        
        min_iter = min_iter or 1
        if max_iter is None or max_iter >= self.data_manager.current_iteration:
            max_iter = self.data_manager.current_iteration - 1 
        
        log.debug('finding trajectory roots for n_iter in [%d,%d]' % (min_iter, max_iter))
        
        roots = []
        for n_iter in xrange(min_iter, max_iter+1):
            # roots always have a p_parent_id less than zero
            p_parent_ids = self.data_manager.get_p_parent_array(n_iter)
            seg_ids = self.data_manager.get_seg_ids_where(n_iter, p_parent_ids < 0)
            segments = self.data_manager.get_segments_by_id(n_iter, seg_ids)
            roots.extend(segments)
        return roots
    
    def get_initial_nodes(self, min_iter = None, max_iter = None, whole_only = False):
        '''Get segments with which to begin a tree walk.  If ``whole_only`` is false (the
        default), return segments which are alive at min_iter or start new trajectories between 
        min_iter and max_iter (inclusive).  If ``whole_only`` is true, return only segments
        which begin new trajectories.'''
        
        min_iter = min_iter or 1
        if max_iter is None or max_iter >= self.data_manager.current_iteration:
            max_iter = self.data_manager.current_iteration - 1
            
        roots = []
        min_iter_roots = set()
        for n_iter in xrange(min_iter, max_iter+1):
            # roots always have a p_parent_id less than zero
            p_parent_ids = self.data_manager.get_p_parent_array(n_iter)
            seg_ids = self.data_manager.get_seg_ids_where(n_iter, p_parent_ids < 0)
            segments = self.data_manager.get_segments_by_id(n_iter, seg_ids)
            roots.extend(segments)
            
            if n_iter == min_iter:
                min_iter_roots.update(seg_ids)
                
        if not whole_only:
            min_iter_segs = self.data_manager.get_segments(min_iter)
            for segment in min_iter_segs:
                if segment.seg_id not in min_iter_roots:
                    roots.append(segment)
        
        return roots
    
    def trace_trajectories(self, min_iter, max_iter, 
                           callable, args=None, kwargs= None, get_state = None, set_state = None,
                           whole_only = False):
        """
        Walk the trajectory tree depth-first, calling
          ``callable(segment, children, history, *args, **kwargs)`` for each segment 
        visited. ``segment`` is the segment being visited, ``children`` is that
        segment's children, ``history`` is the chain of segments leading
        to ``segment`` (not including ``segment``). get_state and set_state are
        used to record and reset, respectively, any state specific to 
        ``callable`` when a new branch is traversed.
        
        If called with whole_only == False (the default), analyze all trajectories 
        alive at min_iter and those newly created between min_iter and max_iter
        (inclusive).  Conversely, if whole_only is true, analyze only trajectories
        which begin between min_iter and max_iter (inclusive), thus analyzing only
        entire trajectories.
        """
        
        args = args or tuple()
        kwargs = kwargs or dict()
        
        self.segments_visited = 0
        
        # Either both or neither of external state getter/setter required
        if (get_state or set_state) and not (get_state and set_state):
            raise ValueError('either both or neither of get_state/set_state must be specified')
            
        # This will grow to contain the maximum trajectory length
        history = numpy.empty((self.history_chunksize,), numpy.object_)
        self.count_segs_in_range(min_iter, max_iter)
        
        #roots = self.get_roots(min_iter,max_iter)
        roots = self.get_initial_nodes(min_iter, max_iter, whole_only)
        #log.debug('%d roots' % len(roots))
                
        for root in roots:
            children = self.data_manager.get_children(root)

            # Visit the root node of each tree unconditionally
            callable(root, children, [], *args, **kwargs)
            self.segments_visited += 1

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
                while node.n_iter < max_iter and len(children):
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
                    children = self.data_manager.get_children(node)

                    # Visit the new node as we descend
                    callable(node, children, history[:len_history], *args, **kwargs)
                    self.segments_visited += 1
                    
        
    def trace_to_root(self, n_iter, seg_id):
        '''Trace the given segment back to its starting point, returning a list of Segment
        objects describing the entire trajectory.'''
        
        segments = []
        segment = self.data_manager.get_segments_by_id(n_iter, [seg_id])[0]
        segments.append(segment)
        while segment.p_parent_id >= 0:
            segment = self.data_manager.get_segments_by_id(segment.n_iter-1, [segment.p_parent_id])[0]
            segments.append(segment)
        return list(reversed(segments))
