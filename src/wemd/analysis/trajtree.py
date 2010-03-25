from sqlalchemy.sql.expression import bindparam
_metaclass__ = type

import numpy

from wemd import Segment, Trajectory
from wemd.data_manager.schema import segmentsTable

import sqlalchemy
from sqlalchemy.sql import select
from sqlalchemy.orm import lazyload

class TrajTree:
    def __init__(self, data_manager, squeeze_data):
        self.data_manager = data_manager
        self.squeeze_data = squeeze_data

        # We only use the data manager for its SQL connection
        # (the ORM is too slow, so we use the more fundamental SQL interface)
        self.bind = data_manager.dbengine.connect()
    
    def _seg_id(self, node):
        try:
            return node.seg_id
        except AttributeError:
            return long(node)
        
    def get_roots(self, max_iter = None):        
        s = select([segmentsTable])
        
        if max_iter:
            s = s.where( (segmentsTable.c.p_parent_id == None)
                          &(segmentsTable.c.n_iter <= max_iter) )
        else:
            s = s.where(segmentsTable.c.p_parent_id == None)
            
        return list(self.bind.execute(s))
    
    def get_children(self, node):
        seg_id = self._seg_id(node)
        
        s = select([segmentsTable]).where(segmentsTable.c.p_parent_id == seg_id)
        return list(self.bind.execute(s))
    
    def _desc_traj(self, node, segs, max_iter, trajs):
        s = select([segmentsTable]).where(segmentsTable.c.p_parent_id == node.seg_id)
        children = self.bind.execute(s).fetchall()
        if not children or node.n_iter == max_iter:
            trajs.append(Trajectory(segs + [node], self.squeeze_data))
        else:
            for child in children:
                self._desc_traj(child, segs + [node], max_iter, trajs)
            
    
    def get_trajectories_r(self, max_iter):
        trajs = []
        for root in self.get_roots(max_iter):
            self._desc_traj(root, [], max_iter, trajs)
        return trajs
         
    def get_trajectories(self, max_iter):
        trajs = []
        
        s = select([segmentsTable]).where((segmentsTable.c.p_parent_id == bindparam('p_parent_id')))
        
        for root in self.get_roots(max_iter):
            state_stack = [{'node': root,
                           'children': list(self.bind.execute(s, {'p_parent_id': root.seg_id})),
                           }]
            
                        
            while state_stack:
                state = state_stack.pop(-1)
                
                # We don't want to descend if we are already at max_iter
                if state['node'].n_iter == max_iter:
                    #trajs.append(state['node'].seg_id)
                    trajs.append([si['node'].seg_id for si in state_stack] + [state['node'].seg_id])
                else:                
                    while state['children']:
                        #state_stack.append(state)
                        node = state['children'].pop(-1)                    
                        children = list(self.bind.execute(s, {'p_parent_id': node.seg_id}))
                        if not children or node.n_iter == max_iter:
                            trajs.append([si['node'].seg_id for si in state_stack] 
                                         + [state['node'].seg_id, node.seg_id])
                            #trajs.append(node.seg_id)
                        else:
                            state_stack.append(state)
                            state = {'node': node,
                                     'children': children}

        return trajs