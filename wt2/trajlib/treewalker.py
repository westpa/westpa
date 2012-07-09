from __future__ import print_function, division; __metaclass__ = type
from itertools import izip
from collections import namedtuple
import numpy
import wemd
from wt2.tool_classes.selected_segs import AllSegmentSelection
from wemd.data_manager import seg_id_dtype

trajnode = namedtuple('trajnode', ('n_iter', 'seg_id'))

def get_roots(start_iter=1, stop_iter = None, segsel = None, data_manager = None):
    '''Return a sequence of (n_iter,seg_id) pairs of segments which initiate new trajectories
    in iterations between start_iter and stop_iter.'''
    
    data_manager = data_manager or wemd.rc.get_data_manager()
    stop_iter = stop_iter or data_manager.current_iteration
    segsel = segsel or AllSegmentSelection(start_iter, stop_iter, data_manager)
    roots = []
    
    # all trajectories in start_iter represent new trajectories for the purposes of analysis
    roots.extend((start_iter, seg_id) for seg_id in segsel.from_iter(start_iter))
    
    # in remaining iterations, only truly newly-initiated trajectories count
    for n_iter in xrange(start_iter+1,stop_iter):
        seg_ids = numpy.array(segsel.from_iter(n_iter))
        parent_ids = numpy.array(data_manager.get_parent_ids(n_iter, seg_ids))
        iter_roots = seg_ids[parent_ids < 0]
        roots.extend((n_iter, seg_id) for seg_id in iter_roots)
        
    return roots

def trace_trajectories(visit, get_visitor_state = None, set_visitor_state = None, cbargs=None, cbkwargs=None,
                       start_iter=1, stop_iter = None, segsel = None,
                       history_chunksize = 1024,
                       data_manager = None):
    '''Trace all trajectories in the weighed ensemble simulation between ``start_iter`` (default 1)
    and ``stop_iter`` (default current iteration), calling 
    ``visit(n_iter, seg_id, *cbargs, **cbkwargs)`` for each segment encountered.
    
    ``get_visitor_state`` and ``set_visitor_state`` are used to manage data for the per-trajectory analysis;
    whenever the tree walk is about to visit one of multiple children, get_visitor_state() is called to save the
    state of the analysis, and whenever the tree walk returns to a previously-visted segment in order to consider its children,
    set_visitor_state(state) is called with the state previously retrieved using get_visitor_state() at that point of the
    calculation. 
    
    If ``segsel`` is provided, it must be a ``SegmentSelection`` which identifies which segments
    to include in the analysis. If absent, all segments between ``start_iter`` and ``stop_iter`` are included.
    
    Returns the number of segments visited.
    '''
    
    if (get_visitor_state or set_visitor_state) and not (get_visitor_state and set_visitor_state):
        raise ValueError('either both or neither of get_visitor_state and set_visitor_state must be specified')
    
    data_manager = data_manager or wemd.rc.get_data_manager()
    
    cbargs = cbargs or ()
    cbkwargs = cbkwargs or {}
    stop_iter = stop_iter or data_manager.current_iteration
    segsel = segsel or AllSegmentSelection(start_iter, stop_iter, data_manager)
    n_visits = 0
    
    roots = get_roots(start_iter, stop_iter, segsel, data_manager)
    print('Examining {:d} roots'.format(len(roots)))
    for (root_n_iter, root_seg_id) in roots:
        state_stack = []
        node = trajnode(root_n_iter, root_seg_id)
        #if node not in segsel:
        #    continue
        
        children = [trajnode(root_n_iter+1, child_id) for child_id in data_manager.get_child_ids(node.n_iter, node.seg_id)]
        state_stack.append({'node': node,
                            'children': children,
                            'xstate': get_visitor_state() if get_visitor_state else None})
        
        while state_stack: 
            state = state_stack.pop(-1)
            node = state['node']
            children = state['children']
            if set_visitor_state:
                set_visitor_state(state['xstate'])
            
            try:
                n_visits += 1
                visit(node.n_iter, node.seg_id, *cbargs, **cbkwargs)
            except StopIteration:
                # Do not descend any further
                del children[:]
                continue
                
            if children:
                node = children.pop(-1)
                #while node not in segsel:
                #    node = children.pop(-1)
                #else:
                    # no children are selected
                #    continue # to next item on the stack
                
                # otherwise, prepare to visit this child
                children = [trajnode(node.n_iter+1, child_id) for child_id in data_manager.get_child_ids(node.n_iter, node.seg_id)]
                state_stack.append({'node': node,
                                    'children': children,
                                    'xstate': get_visitor_state() if get_visitor_state else None})
    return n_visits
