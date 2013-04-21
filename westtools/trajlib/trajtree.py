# Copyright (C) 2013 Matthew C. Zwier
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

__metaclass__ = type
import numpy
import westpa

from itertools import izip
from collections import namedtuple, deque
from west import Segment

TrajID = namedtuple('TrajID', ['n_iter', 'seg_id'])

def all_endpoints(n_iter, iter_group, max_iter = None):
    '''Return all trajectory endpoints (recycled or merged) in the given iteration.
    If the optional ``max_iter`` is supplied, treat all trajectories alive in that iteration
    as endpoints (to facilitate walking the entire tree).'''
    
    assert n_iter > 0    
    if max_iter is not None and n_iter == max_iter:
        return numpy.arange(iter_group['seg_index'].len(), dtype=numpy.int)
    else:
        return (iter_group['seg_index']['endpoint_type'] != Segment.SEG_ENDPOINT_CONTINUES).nonzero()[0]
        

class TrajNode:
    def __init__(self, n_iter, seg_id, parent_id = None, prev = None, next = None, weight=None):
        self.n_iter = int(n_iter)
        self.seg_id = long(seg_id)
        self.parent_id = long(parent_id) if parent_id is not None else None
        self.weight = None
        
        # A reference to the previous segment's node, or None if this is a 
        # trajectory root
        self.prev = prev
        
        # A list of references to child nodes
        self.next = next or set()

    def __hash__(self):
        return id(self)
        
    def __repr__(self):
        return super(TrajNode,self).__repr__()[:-1] + ' n_iter={!r} seg_id={!r} parent_id={!r}>'\
            .format(int(self.n_iter), long(self.seg_id), long(self.parent_id))
            
    def trace(self):
        '''Trace this trajectory back to its root from this node. Returns
        a list of TrajNodes'''
        return list(self.itrace())
    
    def itrace(self):
        '''Trace this trajectory back to its root from this node. Returns
        an iterator over TrajNodes.'''
        node = self
        yield self
        while node.prev:
            yield node.prev
            node = node.prev
            
    def ftrace(self):
        '''Trace this trajectory from its root to this node. Returns a list of TrajNodes'''
        ftrace = self.trace()
        ftrace.reverse()
        return ftrace
    
    def iftrace(self):
        '''Trace this trajectory from its root to this node. Returns an iterator over
        TrajNodes, but must first store the trajectory in forward order.'''
        return iter(self.ftrace())
        
class TrajTree:
    def __init__(self):
        
        # The roots (independent trees) of the trajectories considered here
        # as a mapping from traj_ids to nodes
        self.roots = {}
        
        # The leaves (terminal segments) of the trajectories considered here
        # as a mapping from traj_ids to nodes
        self.leaves = {}
                        
    def __len__(self):
        '''Return the number of segments in this trajectory tree. O(N).'''
        nodes = set(self.leaves.itervalues())
        parents = {node.prev for node in nodes if node.prev is not None}
        nodes.update(parents)
        while parents:
            parents = {node.prev for node in parents if node.prev is not None}
            nodes.update(parents)
        return len(nodes)

    def __iter__(self):
        '''Iterate over nodes, roughly depth-first'''
        state_stack = deque()
        
        for node in self.roots.itervalues():
            yield node
            
            state_stack.append((node, set(node.next)))
            
            while state_stack:
                node, remaining_children = state_stack.pop()
                
                while remaining_children:
                    # Save current state before descending
                    state_stack.append((node, remaining_children))
                    
                    node = remaining_children.pop()
                    remaining_children = set(node.next)
                    
                    # visit new node, then descend on next iteration of this loop
                    yield node
                    
    def write_dot(self, dotfile):
        '''Write a Graphviz Dot input for this tree'''
        
        if isinstance(dotfile, basestring):
            dotfile = open(dotfile, 'wt')
            
        dotfile.write('strict digraph {\n')
        
        for node in self.leaves.itervalues():
            prev = node.prev
            while prev:
                dotfile.write('    "({:d},{:d})" -> "({:d},{:d})";\n'.format(int(prev.n_iter), long(prev.seg_id), 
                                                                         int(node.n_iter), long(node.seg_id)))
                node = prev
                prev = node.prev
                
        dotfile.write('}\n')        
        
def construct_tree(get_matching_segs = all_endpoints, max_iter = None, data_manager = None):
    '''Construct a trajectory tree whose highest iteration is max_iter and whose
    trajectory endpoints in a given iteration are returned by the supplied 
    ``get_matching_segs(n_iter, iter_group)`` function, which takes the iteration 
    number and HDF5 group and returns a sequence of seg_ids.'''
    
    data_manager = data_manager or westpa.rc.get_data_manager()
    max_iter = max_iter or data_manager.current_iteration - 1
        
    tree = TrajTree()
    # a mapping of seg_id in n_iter to parents (seg_id in n_iter-1)
    
    next_iter_nodes = []    
    for n_iter in xrange(max_iter, 0, -1):
        print('building nodes for iteration {:d}'.format(n_iter))
        iter_group = data_manager.get_iter_group(n_iter)
        
        leaf_seg_ids = list(get_matching_segs(n_iter, iter_group))
        leaf_nodes = [TrajNode(n_iter, seg_id) for seg_id in leaf_seg_ids]
        tree.leaves.update({(node.n_iter, node.seg_id): node for node in leaf_nodes})
        
        this_iter_nodes = list(leaf_nodes)
        
        new_nodes = {seg_id: TrajNode(n_iter, seg_id) 
                     for seg_id in set(node.parent_id for node in next_iter_nodes if node.parent_id >= 0)}
        this_iter_nodes.extend(new_nodes.itervalues())
        
        for node in next_iter_nodes:
            if node.parent_id < 0:
                tree.roots[node.n_iter,node.seg_id] = node
            else:
                new_node = new_nodes[node.parent_id]
                node.prev = new_node
                new_node.next.add(node)

        if this_iter_nodes:
            seg_ids = [node.seg_id for node in this_iter_nodes]
            parent_ids = data_manager.get_parent_ids(n_iter, seg_ids)
            weights = data_manager.get_weights(n_iter,seg_ids)
            for node, parent_id, weight in izip(this_iter_nodes, parent_ids, weights):
                node.parent_id = long(parent_id) if parent_id is not None else None
                node.weight = float(weight)
            del seg_ids, weights, parent_ids
            
        next_iter_nodes = this_iter_nodes
        
        del iter_group, leaf_seg_ids, leaf_nodes, this_iter_nodes, new_nodes
        #gc.collect()
        
        
    # roots at the beginning of the simulation are not recorded above, so do that now
    for node in next_iter_nodes:
        assert node.parent_id < 0
        assert node.n_iter == 1
        tree.roots[node.n_iter, node.seg_id] = node
        
    return tree
                
def commonality_matrix(tree):
    '''Calculate a matrix C where C_ij is the number of segments trajectories i and j have
    in common. Returns (leaves, C), where C is the matrix described above and leaves is the
    list of terminal segments which determines the row/column ordering of C (i.e. which 
    trajectory corresponds to which row)'''
    
    n_trajs = len(tree.leaves)
    min_iter = min(root.n_iter for root in tree.roots.itervalues())
    max_iter = max(leaf.n_iter for leaf in tree.leaves.itervalues())
    
    C = numpy.zeros((n_trajs, n_trajs), numpy.min_scalar_type(max_iter-min_iter+1))
    
    leaves = list(tree.leaves.itervalues())
    for i in xrange(n_trajs):
        traj_i = leaves[i].ftrace()
        for j in xrange(i+1):
            traj_j = leaves[j].ftrace()
            
            Cij = 0
            
            for node_i, node_j in izip(traj_i, traj_j):
                if node_i is node_j:
                    Cij += 1
                else:
                    break
            
            C[i,j] = C[j,i] = Cij    

    return leaves, C
                        
          
        
        
        
        
        

        
        
    
    
    