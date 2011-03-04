import logging
log = logging.getLogger(__name__)

import sqlalchemy, sqlalchemy.sql
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, eagerload, lazyload, defer
from sqlalchemy.sql import select, bindparam

COUNT = sqlalchemy.sql.func.count
MIN = sqlalchemy.sql.func.min
MAX = sqlalchemy.sql.func.max

from old_types import OldSegment, WESimIter
import schema

def _n_iter(we_sim_iter):
    try:
        return we_sim_iter.n_iter
    except AttributeError:
        return we_sim_iter
    
class SQLAlchemyDataManager:
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        runtime_config = sim_manager.runtime_config
        runtime_config.require('data.db.url')
    
        # Connect to the database
        self.dbengine = create_engine(runtime_config['data.db.url'],
                                      # don't block on simple queries
                                       isolation_level = 'READ UNCOMMITTED',
                                      # if for some reason we do block, don't die after 5 secs
                                      connect_args={'timeout': 90})
        
        
        # Session factory object
        self.DBSession = sessionmaker(bind=self.dbengine,
                                      autocommit = True,
                                      autoflush = False,
                                      expire_on_commit = False)
        # Transient DB session
        self.dbsession = None
        
        from mappingtable import DictTableIface
        self.meta = DictTableIface(self.dbengine, 
                                   schema.meta_table,
                                   schema.meta_table.c.key_,
                                   schema.meta_table.c.value)

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
                                
    def get_we_sim_iter(self, n_iter):
        self.require_dbsession()
        return self.dbsession.query(WESimIter).filter(WESimIter.n_iter == n_iter).one()
    
    def create_segments(self, we_sim_iter, segments):
        dbsession = self.require_dbsession()
        dbsession.add_all(segments)
                    
    def num_incomplete_segments(self, we_iter):
        self.require_dbsession()
        return self.dbsession.query(OldSegment)\
            .filter(OldSegment.n_iter == _n_iter(we_iter))\
            .filter(OldSegment.status != OldSegment.SEG_STATUS_COMPLETE)\
            .count()

    def get_segment(self, seg_id, **kwargs):
        self.require_dbsession()
        q = self.dbsession.query(OldSegment).filter(OldSegment.seg_id==seg_id)
        return q.options(eagerload(OldSegment.p_parent)).one()
    
    def _get_segment_objects(self, criteria, **kwargs):
        self.require_dbsession()
        q = self.dbsession.query(OldSegment)
        if criteria is not None:
            q = q.filter(criteria)
        
        if kwargs.get('load_p_parent', False):
            q = q.options(eagerload(OldSegment.p_parent))
            
        if kwargs.get('load_parents', False):
            q = q.options(eagerload(OldSegment.parents))
        
        return q.all()
    
    def _get_segment_rows(self, criteria, **kwargs):
        conn = self.dbengine.connect()
        s = select([schema.segments_table])
        if criteria is not None:
            s = s.where(criteria)
            
        try:
            return list(conn.execute(s))
        finally:
            conn.close()
    
    def get_segments(self, 
                     criteria=None, 
                     result_format='objects', # or 'rows' 
                     **kwargs):
        if result_format == 'objects':
            return self._get_segment_objects(criteria, **kwargs)
        elif result_format == 'rows':
            return self._get_segment_rows(criteria, **kwargs)
        else:
            raise ValueError("invalid result format %r; valid choices are 'objects' or 'rows'")
        
    def get_connectivity(self, criteria=None, **kwargs):
        conn = self.dbengine.connect()
        
        s = select([schema.segments_table.c.p_parent_id, schema.segments_table.c.seg_id])
        s = s.where(schema.segments_table.c.p_parent_id != None)
        if criteria is not None:
            s = s.where(criteria)
        try:
            return list(conn.execute(s))
        finally:
            conn.close()
                                
