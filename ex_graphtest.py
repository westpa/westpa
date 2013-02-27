from __future__ import print_function, division
import argparse, numpy
import westpa
from westtools.trajlib import trajtree

parser = argparse.ArgumentParser()
westpa.rc.add_args(parser)
args = parser.parse_args()
westpa.rc.process_args(args)

def get_leaves(n_iter, iter_group):
    seg_index = iter_group['seg_index']
    final_rmsd = iter_group['pcoord'][:,50,1]
    was_recycled = seg_index['endpoint_type'] == 3
    was_bound = final_rmsd < 2.7
    
    return (was_bound & was_recycled).nonzero()[0]

def get_all_leaves(n_iter, iter_group):
    if n_iter == 919:
        return numpy.arange(iter_group['seg_index'].len(), dtype=numpy.int)
    else:
        return (iter_group['seg_index']['endpoint_type'] == 3).nonzero()[0]
    
data_manager = westpa.rc.get_data_manager()
data_manager.open_backing('r')

tree = trajtree.construct_tree(data_manager.current_iteration-1, get_leaves, data_manager)
print('tree has {:d} roots, {:d} leaves, and {:d} segments'.format(len(tree.roots), len(tree.leaves), len(tree)))

all_nodes = list(tree)
print('all_nodes has length {:d}'.format(len(all_nodes)))

leaves_by_branchpoint = {}
leftovers_by_root = {}
for leaf in tree.leaves.itervalues():
    next = None
    for depth, node in enumerate(leaf.itrace()):
        if (depth > 200 and len(node.next) > 1):
            branchpoint_dict = leaves_by_branchpoint.setdefault(node, {})
            leafset = branchpoint_dict.setdefault(next, set())
            leafset.add(leaf)
            #leaves_by_branchpoint.setdefault((node,next), set()).add(leaf)
            break
        elif node.prev is None:
            leftovers_by_root.setdefault(node, set()).add(leaf)
        #    leaves_by_branchpoint.setdefault(node, set()).add(leaf)
        else:
            next = node

#print(leaves_by_branchpoint)
print('there are {:d} branch points leading to approximately independent trajectory bundles'.format(len(leaves_by_branchpoint)))
#for ((branch,node),leaves) in leaves_by_branchpoint.iteritems():
for (branchnode, subtrees) in leaves_by_branchpoint.iteritems():
    print('  branch point {}:{} has {} independent subtree(s)'.format(branchnode.n_iter,branchnode.seg_id, len(subtrees)))
    for (subtree, leaves) in subtrees.iteritems():
        leaves=list(leaves)
        weights = numpy.fromiter((leaf.weight for leaf in leaves), dtype=numpy.float64)
        max_weight_node = leaves[numpy.argmax(weights)]
        print('    subtree {}:{} has {} leaves; max weight ({}) is {}:{}'.format(subtree.n_iter, subtree.seg_id, len(leaves),
                                                                               max_weight_node.weight, max_weight_node.n_iter,
                                                                               max_weight_node.seg_id))
        
        #print('subtree rooted at point ({},{}) contains {:d} leaves'.format(node.n_iter, node.seg_id,len(leaves)))
        #print('  of which ({},{}) has the highest weight ({!r})'.format(max_weight_node.n_iter, max_weight_node.seg_id, 
        #                                                                max_weight_node.weight))
for (root, leaves) in leftovers_by_root.iteritems():
    print('{} trajectories from root {}:{} pruned due to shared history'.format(len(leaves), root.n_iter, root.seg_id))
    
    

#leaves, C = trajtree.commonality_matrix(tree)
#B = numpy.zeros(C.shape, C.dtype)
#F = numpy.zeros(C.shape, numpy.float32)
#for i in xrange(C.shape[0]):
#    F[i,:] = C[i,:] / C[i,i]
#    B[i,:] = C[i,i] - C[i,:]
#rleaves = [(node.n_iter, node.seg_id) for node in leaves]
#pickle.dump(rleaves, open('west_test_common_leaves.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)
#numpy.save('west_test_common.npy', C)
#numpy.save('west_test_common_F.npy', F)
