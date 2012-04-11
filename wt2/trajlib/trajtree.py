__metaclass__ = type
import numpy
import wemd

from itertools import izip
from collections import namedtuple, deque

TrajID = namedtuple('TrajID', ['n_iter', 'seg_id'])

class TrajNode:
    def __init__(self, n_iter, seg_id, parent_id = None, prev = None):
        self.n_iter = n_iter
        self.seg_id = seg_id
        self.parent_id = parent_id
        
        # A reference to the previous segment's node, or None if this is a 
        # trajectory root
        self.prev = prev
        
    def __repr__(self):
        return super(TrajNode,self).__repr__()[:-1] + 'n_iter={!r} seg_id={!r} parent_id={!r}>'\
            .format(int(self.n_iter), long(self.seg_id), self.parent_id )
        

def reverse_traj_iter(terminal_node):
    node = terminal_node
    #print "yielding {!r}".format(node)
    yield node
    while node.prev:
        #print "yielding {!r}".format(node.prev)
        yield node.prev
        node = node.prev
        
def forward_traj_iter(terminal_node):
    traj = list(reverse_traj_iter(terminal_node))
    traj.reverse()
    for node in traj:
        yield node        

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
        n_nodes = 0
        nodes = set(self.leaves.itervalues())
        n_nodes += len(nodes)
        parents = {node.prev for node in nodes if node.prev is not None}
        del nodes
        while parents:
            n_nodes += len(parents)
            parents = {node.prev for node in parents if node.prev is not None}
        return n_nodes
            
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
        
def construct_tree(max_iter, get_matching_segs, data_manager = None):
    '''Construct a trajectory tree whose highest iteration is max_iter and whose
    trajectory endpoints in a given iteration are returned by the supplied 
    ``get_matching_segs(n_iter, iter_group)`` function, which takes the iteration 
    number and HDF5 group and returns a sequence of seg_ids.'''
    
    data_manager = data_manager or wemd.rc.get_data_manager()
        
    tree = TrajTree()
    # a mapping of seg_id in n_iter to parents (seg_id in n_iter-1)
    parent_ids = {}
    # a mapping of seg_id to corresponding nodes
    next_iter_nodes = {}
    
    for n_iter in xrange(max_iter, 0, -1):
        iter_group = data_manager.get_iter_group(n_iter)
        
        this_iter_nodes = {}

        # First, find and record leaves (terminal segments)        
        leaf_seg_ids = list(get_matching_segs(n_iter, iter_group))
        leaf_nodes = [TrajNode(n_iter, seg_id) for seg_id in leaf_seg_ids]
        tree.leaves.update(((node.n_iter, node.seg_id), node) for node in leaf_nodes)
        this_iter_nodes = {node.seg_id: node for node in leaf_nodes}
        
        # Next, create nodes for non-leaf segments in this iteration
        # parent_ids contains a mapping from the seg_ids of the last iteration of this loop (n_iter+1)
        # to parents in this iteration (n_iter)
#        for (child_id, parent_id) in parent_ids.iteritems():
#            next_iter_node = next_iter_nodes[child_id]
#            next_iter_node.parent_id = parent_id
#            
#            if parent_id < 0:
#                tree.roots[n_iter+1,child_id] = next_iter_node
#            else:
#                node = TrajNode(n_iter, parent_id)
#                this_iter_nodes[parent_id] = node
#                next_iter_node.prev = node

        for (seg_id, node) in next_iter_nodes.iteritems():
            parent_id = parent_ids[seg_id]
            node.parent_id = parent_id
            if parent_id < 0:
                tree.roots[node.n_iter,node.seg_id] = node
            else:
                new_node = TrajNode(n_iter, parent_id)
                this_iter_nodes[parent_id] = new_node
                node.prev = new_node

        if this_iter_nodes:                        
            parent_ids = data_manager.get_parent_ids(n_iter, this_iter_nodes.keys())
            assert len(parent_ids) == len(this_iter_nodes)
            #print('n_iter={!r}, parent_ids={!r}'.format(n_iter,parent_ids))
        else:
            parent_ids = {}
        #print('n_iter={!r}, next_iter_nodes={!r}, this_iter_nodes={!r}'.format(n_iter, next_iter_nodes, this_iter_nodes))            
        next_iter_nodes = this_iter_nodes
        
    # roots at the beginning of the simulation are not recorded above, so do that now
    for (child_id, parent_id) in parent_ids.iteritems():
        assert parent_id < 0
        node = next_iter_nodes[child_id]
        assert node.n_iter == 1
        tree.roots[1,child_id] = node
        
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
        traj_i = forward_traj_iter(leaves[i])
        for j in xrange(i):
            traj_j = forward_traj_iter(leaves[j])
            
            for n_iter, (node_i, node_j) in enumerate(izip(traj_i, traj_j)):
                print n_iter, i, j, (node_i.n_iter, node_i.seg_id), (node_j.n_iter, node_j.seg_id),
                if node_i is node_j:
                    print 'match'
                    C[i,j] += 1
                else:
                    print
                    break

    return leaves, C
        
        
          
        
        
        
        
        

        
        
    
    
    