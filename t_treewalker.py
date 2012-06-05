from __future__ import print_function
import sys
from wt2.trajlib.treewalker import trace_trajectories
from wt2.tool_classes.selected_segs import AllSegmentSelection
from wemd.data_manager import WEMDDataManager
from trajtree import TrajTreeSet

data_manager = WEMDDataManager()
data_manager.we_h5filename = 'system.h5'
data_manager.open_backing(mode='r')

segsel = AllSegmentSelection(data_manager = data_manager)
print('There are {:d} segments selected.\n'.format(len(segsel)))

tts = TrajTreeSet(segsel, data_manager)
print('There are {:d} segments in the tree map\n'.format(len(tts)))

visited = set()
def visit(n_iter, seg_id, visited):
    raise StopIteration()
    if (n_iter,seg_id) in visited:
        raise ValueError('{:d}:{:d} already seen'.format(n_iter,seg_id))
    else:
        visited.add((n_iter,seg_id))
        print('\rVisited {:>6d}:{:<6d}'.format(n_iter,seg_id),end='')
        if len(visited) % 10:
            sys.stdout.flush()
    

#n_visited = trace_trajectories(visit, cbargs=(visited,), segsel=segsel, data_manager=data_manager)
#print('\nVisited {:d} segments.'.format(n_visited))
    




