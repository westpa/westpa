import schema
import sqlalchemy
from sqlalchemy import select
from sqlalchemy.exceptions import OperationalError

import cPickle as pickle

import logging
log = logging.getLogger(__name__)

db_version_key = 'db.schema_version'

_updaters = []

def get_schema_version(engine):
    try:
        rsl = engine.execute("SELECT key_, value FROM meta WHERE key_='%s'"
                             %  db_version_key)
    except OperationalError, e:
        if 'table' in str(e).lower():
            return 0
        else:
            raise
    else:
        try:
            row = rsl.fetchone()
            if row is None: return 0
            return pickle.loads(str(row[1]))
        finally:
            rsl.close()
                                
def update_schema(engine):
    global schema
    metadata = schema.metadata
    
    for updater_class in _updaters:
        updater = updater_class()
        if updater.needs_update(schema, metadata, engine):
            updater.apply_update(schema, metadata, engine)

class SchemaUpdate(object):
    target_version = None
    
    def __init__(self):
        self.log = logging.getLogger(self.__class__.__name__)
        
    def create_tables(self, schema, metadata, engine, tables):
        self.log.info('creating new tables')
        metadata.create_all(bind=engine, tables=tables, checkfirst=True)
        
    def drop_tables(self, schema, metadata, engine, tables):
        self.log.info('dropping outdated tables')
        metadata.drop_all(bind=engine, tables=tables, checkfirst=True)
    
    def needs_update(self, schema, metadata, engine):
        assert self.target_version is not None
        if get_schema_version(engine) < self.target_version:
            return True
    
    def apply_update(self, schema, metadata, engine):
        self.pre_update(schema, metadata, engine)
        self.update_schema(schema, metadata, engine)
        self.post_update(schema, metadata, engine)
        self.log.info('database schema updated to version %s' % self.target_version)

    def pre_update(self, schema, metadata, engine):
        pass
    
    def update_schema(self, schema, metadata, engine):
        raise NotImplementedError
    
    def post_update(self, schema, metadata, engine):
        pass
    
class Update_0_1(SchemaUpdate):
    target_version = 1
    
    def update_schema(self, schema, metadata, engine):
        new_tables = ['meta', 'traj_tree']
        
        conn = engine.connect()
        self.create_tables(schema, metadata, engine, new_tables)
                
        sel_iter_0_segs = select([schema.segments_table.c.seg_id],
                                 schema.segments_table.c.n_iter == 0) 
        
        trans = conn.begin()
        try:
            log.info('removing references to iteration 0 segments')
            engine.execute(schema.segment_lineage_table.delete(schema.segment_lineage_table.c.seg_id.in_(sel_iter_0_segs)))
            log.info('removing iteration 0 segments')
            engine.execute(schema.segments_table.delete(schema.segments_table.c.n_iter == 0))
            log.info('setting iteration 1 parents to NULL')
            engine.execute(schema.segments_table.update(schema.segments_table.c.n_iter == 1,
                                                       {'p_parent_id': None}))
            trans.commit()
        except:
            trans.rollback()
            raise
        
        conn.close()

    def post_update(self, schema, metadata, engine):
        metaInsert = schema.meta_table.insert()
        engine.execute(metaInsert, {'key_': db_version_key, 
                                    'value': self.target_version})

class Update_1_2(SchemaUpdate):
    target_version = 2
    
    def update_schema(self, schema, metadata, engine):
        self.drop_tables(schema, metadata, engine, ['traj_tree'])
        
    def post_update(self, schema, metadata, engine):
        self.log.warning('endpoint data not available for this simulation')
    
_updaters.append(Update_0_1)
_updaters.append(Update_1_2)
