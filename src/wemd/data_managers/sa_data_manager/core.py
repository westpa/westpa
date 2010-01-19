from wemd.data_managers import DataManagerBase
from wemd.core import Segment, WESimIter
import schema
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, eagerload, lazyload, defer
from copy import copy

import logging
log = logging.getLogger(__name__)

def _n_iter(we_sim_iter):
    try:
        return we_sim_iter.n_iter
    except AttributeError:
        return we_sim_iter

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
        self.close_dbsession()

    def create_we_sim_iter(self, we_sim_iter):
        self.require_dbsession()
        self.dbsession.add(we_sim_iter)
        self.dbsession.flush()        
    
    def update_we_sim_iter(self, we_sim_iter):
        self.require_dbsession()
        # Ensure an update of the data field.  Since this is a shallow copy,
        # this shouldn't be a memory explosion even with a lot of data
        we_sim_iter.addtl_data = copy(we_sim_iter.data)
        self.dbsession.add(we_sim_iter)
        self.dbsession.flush()
    
    def get_we_sim_iter(self, n_iter):
        self.require_dbsession()
        return self.dbsession.query(WESimIter).filter(WESimIter.n_iter == n_iter).one()
    
    def create_segments(self, we_sim_iter, segments):
        self.require_dbsession()
        dbsession = self.dbsession
        
        dbsession.begin()
        try:
            dbsession.add_all(segments)
        except:
            dbsession.rollback()
            raise
        else:
            dbsession.commit()
    
    def update_segments(self, we_sim_iter, segments, **kwargs):
        self.require_dbsession()
        dbsession = self.dbsession
        
        dbsession.begin()
        try:
            for segment in segments:
                segment.pcoord = copy(segment.pcoord)
                segment.data = copy(segment.data)
                try:
                    segment.p_parent.pcoord = copy(segment.p_parent.pcoord)
                    segment.p_parent.data = copy(segment.p_parent.data)
                except AttributeError:
                    pass
                dbsession.merge(segment)
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
        n_iter = _n_iter(we_iter)
        return self.dbsession.query(Segment)\
            .filter( (Segment.n_iter == n_iter)
                    &(Segment.seg_id == seg_id))\
            .options(eagerload(Segment.p_parent)).one()
    
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
    