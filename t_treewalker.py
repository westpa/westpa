from __future__ import print_function
import sys
from wt2.trajlib.treewalker import trace_trajectories
from wt2.tool_classes.selected_segs import AllSegmentSelection
from west.data_manager import WESTDataManager
from trajtree import TrajTreeSet
from trajtree.trajtree import FakeTrajTreeSet

data_manager = WESTDataManager()
data_manager.we_h5filename = 'system.h5'
data_manager.open_backing(mode='r')

segsel = AllSegmentSelection(data_manager = data_manager)
print('There are {:d} segments selected.\n'.format(len(segsel)))

#tts = TrajTreeSet(segsel, data_manager)
tts = FakeTrajTreeSet()
print('There are {:d} segments in the tree map'.format(len(tts)))
print('  Roots: {:d}, leaves: {:d}'.format(tts.count_roots(), tts.count_leaves()))

visited = set()
def visit(n_iter, seg_id, weight, visited):
    #raise StopIteration()
    #print('Visiting {:>6d}:{:<6d}'.format(n_iter,seg_id))
    if seg_id == 9: raise StopIteration
    else: print('visiting {:d}:{:d}'.format(n_iter,seg_id))
    if (n_iter,seg_id) in visited:
        raise ValueError('{:d}:{:d} already seen'.format(n_iter,seg_id))
    else:
        visited.add((n_iter,seg_id))
    
n_visited = tts.trace_trajectories(visit, vargs=(visited,))
print('Visited {:d} segments'.format(n_visited))
#n_visited = trace_trajectories(visit, cbargs=(visited,), segsel=segsel, data_manager=data_manager)
#print('\nVisited {:d} segments.'.format(n_visited))
    




