import time
from collections import deque

import numpy

import schema
import sqlalchemy
from sqlalchemy.sql import select, bindparam, and_
from sqlalchemy.sql.functions import count as COUNT, max as MAX, min as MIN, sum as SUM

import logging
log = logging.getLogger(__name__)

class ResultSetIterator(object):
    def __init__(self, rslset):
        self.rslset = rslset
        
    def __iter__(self):
        return self
    
    def next(self):
        row = self.rslset.fetchone()
        if row is None:
            raise StopIteration
        else:
            return row

class TrajectoryPath(object):
    def __init__(self, seg_ids, weight, pcoord, 
                 cputime = None, walltime = None,
                 startdate = None, enddate = None,
                 data_first = None, data_last = None):
        self.seg_ids = seg_ids
        self.weight = weight
        self.pcoord = pcoord
        self.cputime = cputime
        self.walltime = walltime
        self.startdate = startdate
        self.enddate = enddate
        self.data_first = data_first
        self.data_last = data_last
        
class SegNode(object):
    __slots__ = ('seg_id',
                 'lt', 'rt')
    def __init__(self, seg_id,
                 lt = None, rt=None):
        self.seg_id = seg_id
        self.lt = lt
        self.rt = rt

class TrajectoryMap(object):  
    TT_INFO_KEY = 'traj_tree_info'  
    def __init__(self, bind, meta, squeeze_data = True):
        self.bind = bind
        self.meta = meta
        self.squeeze_data = squeeze_data
            
    def _find_roots(self, min_iter, max_iter):
        """Query the database and return an iterable of `SegNode` objects 
        representing roots of trajectory trees (i.e. those that have no 
        parents)"""
        
        root_sel = select([schema.segmentsTable.c.seg_id],
                            (schema.segmentsTable.c.p_parent_id == None)
                           &(schema.segmentsTable.c.n_iter >= min_iter)
                           &(schema.segmentsTable.c.n_iter <= max_iter))
        
        
        return iter(SegNode(row.seg_id) 
                    for row in ResultSetIterator(self.bind.execute(root_sel)))
        
    def _get_max_iter(self):
        try:
            info = self.meta[self.TT_INFO_KEY]
        except KeyError:
            return None
        else:
            return info['max_iter']
        
    def _get_min_iter(self):
        try:
            info = self.meta[self.TT_INFO_KEY]
        except KeyError:
            return None
        else:
            return info['min_iter']
        
    def _get_min_max_iter(self):
        try:
            info = self.meta[self.TT_INFO_KEY]
        except KeyError:
            return (None, None)
        else:
            return (info['min_iter'], info['max_iter'])
        
    def is_updated(self, max_iter, min_iter = 1):
        try:
            info = self.meta[self.TT_INFO_KEY]
        except KeyError:
            return False
        else:
            return (info['max_iter'] == max_iter and info['min_iter'] == min_iter)
            
        
    max_iter = property(_get_max_iter, None, None, None)
    min_iter = property(_get_min_iter, None, None, None)
            
    def clear(self):
        log.info('clearing trajectory tree')
        self.bind.execute(schema.trajTreeTable.delete())
        try:
            del self.meta[self.TT_INFO_KEY]
        except KeyError:
            pass
        
    def update(self, max_iter, min_iter=1):
        self.clear()
        
        sel_by_parent = select([schema.segmentsTable.c.seg_id],
                               (schema.segmentsTable.c.p_parent_id == bindparam('p_parent_id'))
                               &(schema.segmentsTable.c.n_iter >= min_iter)
                               &(schema.segmentsTable.c.n_iter <= max_iter))
        
        tt_insert = schema.trajTreeTable.insert()
        nroots = self._get_nroots(max_iter, min_iter)
        
        label = 0
        for (iroot, root) in enumerate(self._find_roots(min_iter, max_iter)):
            log.info('mapping root %d of %d' % (iroot+1, nroots))
            
            node_stack = [root]
            children_stack = [[SegNode(row.seg_id) for row in
                               self.bind.execute(sel_by_parent,
                                                 {'p_parent_id': root.seg_id}).fetchall()]]
            
            while node_stack:
                node = node_stack.pop(-1)
                children = children_stack.pop(-1)
                
                while children:
                    node_stack.append(node)
                    if node.lt is None:
                        label += 1
                        node.lt = label
                    
                    node = children.pop(0)
                    children_stack.append(children)
                    children = [SegNode(row.seg_id) for row in
                                self.bind.execute(sel_by_parent, {'p_parent_id': node.seg_id}).fetchall()]
                else:
                    # For leaf nodes
                    if node.lt is None:
                        label += 1
                        node.lt = label

                assert node.lt is not None
                label += 1
                node.rt = label
                self.bind.execute(tt_insert, {'seg_id': node.seg_id,
                                              'lt': node.lt,
                                              'rt': node.rt})
        
        self.meta[self.TT_INFO_KEY] = {'min_iter': min_iter,
                                       'max_iter': max_iter}
    
    def check_update(self, max_iter, min_iter=1):
        """Update the trajectory tree if not already up-to-date."""
        try:
            tree_info = self.meta[self.TT_INFO_KEY]
        except KeyError:
            log.info('no trajectory tree found; building')
            self.update(max_iter, min_iter)
        else:
            if tree_info['max_iter'] == max_iter and tree_info['min_iter'] == min_iter:
                log.info('trajectory tree is up-to-date')
                return
            else:
                log.info('trajectory tree is out-of-date; rebuilding')
                self.update(max_iter, min_iter)

    def _get_roots(self):
        min_iter, max_iter = self._get_min_max_iter()
        return self._find_roots(min_iter, max_iter)
    
    def _get_nroots(self, max_iter=None, min_iter=1):
        if max_iter is None:
            min_iter, max_iter = self._get_min_max_iter()
        nroot_sel = select([COUNT(schema.segmentsTable.c.seg_id)],
                            (schema.segmentsTable.c.p_parent_id == None)
                           &(schema.segmentsTable.c.n_iter >= min_iter)
                           &(schema.segmentsTable.c.n_iter <= max_iter))
        return self.bind.execute(nroot_sel).scalar()

    
    def _get_leaves(self):
        sel_leaves = select([schema.trajTreeTable.c.seg_id],
                            schema.trajTreeTable.c.rt == schema.trajTreeTable.c.lt + 1)
        return iter(SegNode(row.seg_id) for row in 
                    ResultSetIterator(self.bind.execute(sel_leaves)))
        
    def _get_nleaves(self):
        sel_nleaves = select([COUNT(schema.trajTreeTable.c.seg_id)],
                            schema.trajTreeTable.c.rt == schema.trajTreeTable.c.lt + 1)
        return self.bind.execute(sel_nleaves).scalar()
    
    nroots = property(_get_nroots, None, None, 'Return the number of tree roots')
    nleaves = property(_get_nleaves, None, None, 'Return the number of tree leaves')
    ntrajs = nleaves
    roots = property(_get_roots, None, None, 'Return an iterator over tree roots')
    leaves = property(_get_leaves, None, None, 'Return an iterator over tree leaves')
    
    def get_trajectory(self, leaf_id):
        try:
            leaf_id = leaf_id.seg_id
        except AttributeError:
            leaf_id = long(leaf_id)
        
        parent_tbl = schema.trajTreeTable.alias('parent')
        node_tbl = schema.trajTreeTable.alias('node')
                
        sel_path = select([parent_tbl.c.seg_id],
                          and_(node_tbl.c.lt.between(parent_tbl.c.lt, parent_tbl.c.rt), 
                               node_tbl.c.seg_id == leaf_id),
                          from_obj=[parent_tbl, node_tbl])\
                   .order_by(parent_tbl.c.lt)

#        sel_aggr = select([COUNT(schema.segmentsTable.c.seg_id),
#                           SUM(schema.segmentsTable.c.cputime),
#                           SUM(schema.segmentsTable.c.walltime),
#                           MIN(schema.segmentsTable.c.startdate),
#                           MAX(schema.segmentsTable.c.enddate)],
#                           schema.segmentsTable.c.seg_id.in_(sel_path),
#                           from_obj = [schema.segmentsTable])
        
        sel_details = select([schema.segmentsTable],
                             schema.segmentsTable.c.seg_id.in_(sel_path))
        
        q_start_time = time.clock()
        segments = self.bind.execute(sel_details).fetchall()
        n_segs = len(segments)
        q_end_time = time.clock()
        log.debug('retrieved %d-segment path in %g seconds'
                  % (n_segs, q_end_time - q_start_time))
        
        # Determine shape of progress coordinate and weight arrays
        model_segment = segments[0]
        n_md_steps = model_segment.pcoord.shape[0]
        ndim_pcoord = model_segment.pcoord.ndim - 1
        
        if self.squeeze_data:
            pcoord_shape = ((n_md_steps-1)*(n_segs-1) + n_md_steps,) \
                            + model_segment.pcoord.shape[1:]
        else:
            pcoord_shape = (n_md_steps * n_segs,) \
                            + model_segment.pcoord.shape[1:]

        weight_shape = pcoord_shape[0:1]
        
        seg_ids = numpy.empty((n_segs,), numpy.uint64)
        pcoord = numpy.empty(pcoord_shape, model_segment.pcoord.dtype)
        weight = numpy.empty(weight_shape, numpy.float64)
        
        seg_ids[0] = segments[0].seg_id
        weight[0:n_md_steps] = segments[0].weight
        pcoord[0:n_md_steps] = segments[0].pcoord
        
        for (oiseg, seg) in enumerate(segments[1:]):
            iseg = oiseg+1
            seg_ids[iseg] = seg.seg_id
            
            if self.squeeze_data:
                lb = (n_md_steps-1)*(iseg-1) + n_md_steps
                ub = (n_md_steps-1)*iseg     + n_md_steps
                pco = 1
                assert (seg.pcoord[0] == pcoord[lb-1]).all()
            else:
                lb = (n_md_steps)*(iseg-1)
                ub = (n_md_steps)*iseg
                pco = 0
            
            weight[lb:ub] = seg.weight
            pcoord[lb:ub] = seg.pcoord[pco:]
        
        traj = TrajectoryPath(seg_ids, 
                              weight = weight,
                              pcoord = pcoord,
                              data_first = segments[0].data,
                              data_last = segments[-1].data
                              )
        return traj
        
    def itertrajs(self):
        return iter(self.get_trajectory(leaf_node.seg_id)
                    for leaf_node in self.leaves)
    
    trajs = property(itertrajs, None, None, 'An iterator over trajectory paths')
    
    
