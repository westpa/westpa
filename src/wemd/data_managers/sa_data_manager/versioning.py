import schema
import sqlalchemy
from sqlalchemy import select
from sqlalchemy.exceptions import OperationalError

import cPickle as pickle

import logging
log = logging.getLogger(__name__)

_db_version_key = 'db.schema_version'

_updaters = []

def get_schema_version(engine):
    try:
        rsl = engine.execute("SELECT key_, value FROM meta WHERE key_='%s'"
                             %  _db_version_key)
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
        log = logging.getLogger(self.__class__.__name__)
    
    def needs_update(self, schema, metadata, engine):
        assert self.target_version is not None
        if get_schema_version(engine) < self.target_version:
            return True
    
    def apply_update(self, schema, metadata, engine):
        self.pre_update(schema, metadata, engine)
        self.update_schema(schema, metadata, engine)
        self.post_update(schema, metadata, engine)
        log.info('database schema updated to version %s' % self.target_version)

    def pre_update(self, schema, metadata, engine):
        pass
    
    def update_schema(self, schema, metadata, engine):
        raise NotImplementedError
    
    def post_update(self, schema, metadata, engine):
        pass
    
class Update_0_1(SchemaUpdate):
    target_version = 1
    
    def update_schema(self, schema, metadata, engine):
        new_tables = [schema.metaTable, schema.trajTreeTable]
        metadata.create_all(bind=engine, tables=new_tables, checkfirst=True)

    def post_update(self, schema, metadata, engine):
        metaInsert = schema.metaTable.insert()
        engine.execute(metaInsert, {'key_': _db_version_key, 
                                    'value': self.target_version})

_updaters.append(Update_0_1)
