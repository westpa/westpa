from __future__ import print_function, division; __metaclass__=type
from itertools import imap
import numpy
from wemd import Segment

import logging
log = logging.getLogger(__name__)

class TrajTree:
    def __init__(self, data_manager):
        self.data_manager = data_manager
        
        self.history_chunksize = 32
        
        self.current_iteration = data_manager.current_iteration
        self._iter_groups = {n_iter: data_manager.get_iter_group(n_iter)
                             for n_iter in xrange(1,self.current_iteration+1)}
        self._seg_id_ranges = {n_iter: numpy.arange(0,self._iter_groups[n_iter]['seg_index'].shape[0],
                                                    dtype = numpy.uintp) 
                               for n_iter in xrange(1,self.current_iteration+1)}
        
        # Mapping of iteration -> data from that iteration
        self._seg_indices = dict()
        self._parent_arrays = dict()
        self._p_parent_arrays = dict()
        self._pcoord_arrays = dict()
        
        self.segments_to_visit = 0
        self.segments_visited = 0
        
    def _get_seg_index(self, n_iter):
        try:
            return self._seg_indices[n_iter]
        except KeyError:
            seg_index = self._seg_indices[n_iter] = self._iter_groups[n_iter]['seg_index'][...]
            return seg_index
        
    def _get_parent_array(self, n_iter):
        try:
            return self._parent_arrays[n_iter]
        except KeyError:
            parent_array = self._parent_arrays[n_iter] = self._iter_groups[n_iter]['parents'][...]
            return parent_array
                
    def _get_p_parent_array(self, n_iter):
        try:
            return self._p_parent_arrays[n_iter]
        except KeyError:
            parent_offsets = self._get_seg_index(n_iter)['parents_offset'][:]
            p_parent_array = self._get_parent_array(n_iter)[parent_offsets]
            self._p_parent_arrays[n_iter] = p_parent_array
            return p_parent_array
                
    def _get_pcoord_array(self, n_iter):
        try:
            return self._pcoord_arrays[n_iter]
        except KeyError:
            pcoords = self._pcoord_arrays[n_iter] = self._iter_groups[n_iter]['pcoord'][...]
            return pcoords
        
    def _get_segments_by_id(self, n_iter, seg_ids):
        '''Get segments from the data manager, employing caching'''
        if len(seg_ids) == 0: return []
        
        seg_index  = self._get_seg_index(n_iter)
        all_parent_ids = self._get_parent_array(n_iter)
        
        # Believe it or not, because we're walking over a tree with shared history, and 
        # the segment pcoords are partially shared, it's more memory efficient to 
        # load the entire pcoord array into memory for a given iteration, and let
        # numpy shared views on that array populate the pcoord fields
        pcoords = self._get_pcoord_array(n_iter)
        
        segments = []
        # seg_ids must be a list, or else the pcoords_by_seg dereference will fail
        seg_ids = list(seg_ids)
        
        for seg_id in seg_ids:
            row = seg_index[seg_id]
            parents_offset = row['parents_offset']
            n_parents = row['n_parents']
            segment = Segment(seg_id = seg_id,
                              n_iter = n_iter,
                              status = row['status'],
                              n_parents = n_parents,
                              endpoint_type = row['endpoint_type'],
                              walltime = row['walltime'],
                              cputime = row['cputime'],
                              weight = row['weight'])
            segment.pcoord = pcoords[seg_id,...]
            parent_ids = all_parent_ids[parents_offset:parents_offset+n_parents]
            segment.p_parent_id = long(parent_ids[0])
            segment.parent_ids = set(imap(long,parent_ids))
            segments.append(segment)
        
        return segments
        
        
    def _get_children(self, segment):        
        # get_children gets called for every segment being processed, always,
        # which means that the segment table and parent table get pulled from
        # HDF5 TWICE per segment if data_manager.get_children() is called
        p_parents = self._get_p_parent_array(segment.n_iter+1)
        seg_ids = self._seg_ids_where(segment.n_iter+1, p_parents == segment.seg_id)
        return self._get_segments_by_id(segment.n_iter+1, seg_ids)
    
        
    def _seg_ids_where(self, n_iter, bool_array):
        #seg_index = self._get_seg_index(n_iter)
        #all_seg_ids = numpy.arange(len(seg_index), dtype=numpy.uintp)
        #seg_ids = all_seg_ids[bool_array]
        seg_ids = self._seg_id_ranges[n_iter][bool_array]
        try:
            len(seg_ids)
        except TypeError:
            return [seg_ids]
        else:
            return seg_ids
        
    def count_segs_in_range(self, min_iter = None, max_iter = None):
        min_iter = min_iter or 1
        if max_iter is None or max_iter >= self.data_manager.current_iteration:
            max_iter = self.data_manager.current_iteration - 1 
        
        self.segments_to_visit = 0
        for n_iter in xrange(min_iter, max_iter+1):
            seg_index = self._get_seg_index(n_iter)
            self.segments_to_visit += len(seg_index)        
                
    def get_roots(self, min_iter = None, max_iter = None):
        '''Get segments which start new trajectories.  If min_iter or max_iter is specified, restrict the
        set of iterations within which the search is conducted.'''
        
        min_iter = min_iter or 1
        if max_iter is None or max_iter >= self.data_manager.current_iteration:
            max_iter = self.data_manager.current_iteration - 1 
        
        roots = []
        for n_iter in xrange(min_iter, max_iter+1):
            # roots always have a p_parent_id less than zero
            p_parent_ids = self._get_p_parent_array(n_iter)
            seg_ids = self._seg_ids_where(n_iter, p_parent_ids < 0)
            segments = self._get_segments_by_id(n_iter, seg_ids)
            #segments = self.data_manager.get_segments(n_iter)
            #segments = [segment for segment in segments if segment.p_parent_id < 0]
            roots.extend(segments)
        return roots
    
    def trace_trajectories(self, min_iter, max_iter, 
                           callable, args=None, kwargs= None, get_state = None, set_state = None):
        """
        Walk the trajectory tree depth-first, calling
          ``callable(segment, children, history, *args, **kwargs)`` for each segment 
        visited. ``segment`` is the segment being visited, ``children`` is that
        segment's children, ``history`` is the chain of segments leading
        to ``segment`` (not including ``segment``). get_state and set_state are
        used to record and reset, respectively, any state specific to 
        callable when a new branch is traversed.
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
        
        for root in self.get_roots(min_iter, max_iter):
            children = self._get_children(root)

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
                    children = self._get_children(node)

                    # Visit the new node as we descend
                    callable(node, children, history[:len_history], *args, **kwargs)
                    self.segments_visited += 1
                    
