from copy import copy
import logging
log = logging.getLogger(__name__)

import numpy

import sqlalchemy, sqlalchemy.sql
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, eagerload, lazyload, defer
from sqlalchemy.sql import select, bindparam
COUNT = sqlalchemy.sql.func.count
MIN = sqlalchemy.sql.func.min
MAX = sqlalchemy.sql.func.max

from wemd.core import Segment, WESimIter, Trajectory
import schema, versioning

def _n_iter(we_sim_iter):
    try:
        return we_sim_iter.n_iter
    except AttributeError:
        return we_sim_iter
    
class DataManagerBase(object):
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config

class SQLAlchemyDataManager(DataManagerBase):
    def __init__(self, runtime_config):
        super(SQLAlchemyDataManager,self).__init__(runtime_config)
        runtime_config.require('data.db.url')
    
        # Connect to the database
        self.dbengine = create_engine(runtime_config['data.db.url'])
        
        # Session factory object
        self.DBSession = sessionmaker(bind=self.dbengine,
                                      autocommit = True,
                                      autoflush = False,
                                      expire_on_commit = False)
        # Transient DB session
        self.dbsession = None
        
        from mappingtable import DictTableIface
        self.meta = DictTableIface(self.dbengine, schema.metaTable,
                                   schema.metaTable.c.key_,
                                   schema.metaTable.c.value)

    def new_dbsession(self):
        assert self.dbsession is None
        self.dbsession = self.DBSession()
        
    def require_dbsession(self):
        if self.dbsession is None:
            self.dbsession = self.DBSession()
        return self.dbsession
        
    def close_dbsession(self):
        self.dbsession.close()
        self.dbsession = None
                        
    def prepare_backing(self, sim_config):
        self.new_dbsession()
        log.info('creating database tables')
        schema.metadata.create_all(bind = self.dbengine)
        self.meta[versioning.db_version_key] = schema.version
        self.close_dbsession()

    def create_we_sim_iter(self, we_sim_iter):
        self.require_dbsession()
        self.dbsession.add(we_sim_iter)
        self.dbsession.flush()        
    
    def update_we_sim_iter(self, we_sim_iter):
        self.require_dbsession()
        # Ensure an update of the data field.  Since this is a shallow copy,
        # this shouldn't be a memory explosion even with a lot of data
        # All the gymnastics are because the data array can contain numpy
        # arrays, which throw exceptions when compared using ==/!= (used
        # during the merge() call by SQLAlchemy)
        data = copy(we_sim_iter.data)
        we_sim_iter.data = None
        we_sim_iter = self.dbsession.merge(we_sim_iter)
        self.dbsession.add(we_sim_iter)
        self.dbsession.flush()
        we_sim_iter.data = data
        self.dbsession.flush()
    
    def get_we_sim_iter(self, n_iter):
        self.require_dbsession()
        return self.dbsession.query(WESimIter).filter(WESimIter.n_iter == n_iter).one()
    
    def create_segments(self, we_sim_iter, segments):
        dbsession = self.require_dbsession()
        dbsession.add_all(segments)
    
    def update_segments(self, we_sim_iter, segments, **kwargs):
        dbsession = self.require_dbsession()
        
        dbsession.begin(subtransactions=True)
        try:
            for segment in segments:
                pcoord = copy(segment.pcoord)
                data = copy(segment.data)
                segment.pcoord = None
                segment.data = None
                segment = dbsession.merge(segment)
                dbsession.add(segment)
                dbsession.flush()
                segment.pcoord = pcoord
                segment.data = data
                dbsession.flush()
        except:
            dbsession.rollback()
            raise
        else:
            dbsession.commit()
                
    def num_incomplete_segments(self, we_iter):
        self.require_dbsession()
        return self.dbsession.query(Segment)\
            .filter(Segment.n_iter == _n_iter(we_iter))\
            .filter(Segment.status != Segment.SEG_STATUS_COMPLETE)\
            .count()

    def get_segment(self, we_iter, seg_id, **kwargs):
        self.require_dbsession()
        q = self.dbsession.query(Segment).filter(Segment.seg_id==seg_id)
        if we_iter is not None:
            q = q.filter(Segment.n_iter == _n_iter(we_iter))
        return q.options(eagerload(Segment.p_parent)).one()
    
    def get_segments(self, we_iter, **kwargs):
        self.require_dbsession()
        n_iter = _n_iter(we_iter)
        return self.dbsession.query(Segment)\
            .filter(Segment.n_iter == n_iter)\
            .options(eagerload(Segment.p_parent)).all()
    
    def get_prepared_segments(self, we_iter, **kwargs):
        self.require_dbsession()
        n_iter = _n_iter(we_iter)
        return self.dbsession.query(Segment) \
            .filter( (Segment.n_iter == n_iter)
                    &(Segment.status == Segment.SEG_STATUS_PREPARED))\
            .options(eagerload(Segment.p_parent)).all()
            
    def get_schema(self):
        return schema
    
    def get_schema_version(self):
        return versioning.get_schema_version(self.dbengine)
    
    def get_trajectory(self, leaf_id, squeeze_data = True):
        self.require_dbsession()
        
        seg = self.get_segment(None, leaf_id)
        segs = [seg]
        while seg.p_parent:
            seg = seg.p_parent
            segs.append(seg)
            
        traj = Trajectory(segs, squeeze_data)
        return traj

    def get_trajectories(self, leaf_ids, squeeze_data = True):
        self.require_dbsession()
        seg_id_sel = select([Segment.seg_id, Segment.n_iter])
        seg_sel = select([schema.segmentsTable])
                          
        mit_sel = select([MAX(Segment.n_iter)],
                         Segment.seg_id.in_(leaf_ids))
        max_iter = self.dbsession.execute(mit_sel).fetchone()[0]
        
        traj_data_shape = (max_iter, len(leaf_ids))
        ncols = len(leaf_ids)
        seg_ids = numpy.zeros(traj_data_shape, numpy.uint64)
        segments = numpy.empty(traj_data_shape, numpy.object_)
        segments[...] = None
        
        traj_cols = dict((leaf_id, i) for (i,leaf_id) in enumerate(leaf_ids))
        
        # Assign endpoints
        rsl = self.dbsession.execute(seg_id_sel.where(Segment.seg_id.in_(leaf_ids)))
        for dbrow in rsl:
            trajrow = dbrow.n_iter-1
            trajcol = traj_cols[dbrow.seg_id]
            seg_ids[trajrow, trajcol] = dbrow.seg_id
        
        for iiter in xrange(max_iter-1, -1, -1):
            # What segments exist at this iteration?
            row_ids = seg_ids[iiter,:]
            active_ids = set(long(row_id) for row_id in row_ids[row_ids != 0])
            
            rsl = self.dbsession.execute(seg_sel.where(Segment.seg_id.in_(active_ids)))
            for dbrow in rsl:
                for icol in xrange(0, ncols): 
                    if seg_ids[iiter,icol] == dbrow.seg_id:
                        segments[iiter, icol] = dbrow
                        if iiter > 0:
                            seg_ids[iiter-1, icol] = dbrow.p_parent_id or 0
            rsl.close()
        
        # We now have a matrix whose columns are trajectories
        trajs = []
        for icol in xrange(0, ncols):
            trajcol = segments[:,icol]
            trajsegs = [seg for seg in trajcol if seg is not None]
            trajs.append(Trajectory(trajsegs, squeeze_data))
            
        return trajs
        
        
        
        
        
        
    
