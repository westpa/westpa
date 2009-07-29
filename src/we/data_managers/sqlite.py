import os, errno
try:
    import cPickle as pickle
except ImportError:
    import pickle

try:
    import pysqlite as sqlitedb    
except ImportError:
    import sqlite3 as sqlitedb

import numpy
    
from we.data_managers import DataManagerBase
from we.core.segments import Segment

import logging
log = logging.getLogger('we.data_managers.sqlite')

LOCK_TIMEOUT = 120

def pickle_to_blob(obj):
    return buffer(pickle.dumps(obj, pickle.HIGHEST_PROTOCOL))

def unpickle_from_blob(field):
    if field is None: 
        return None
    else:
        return pickle.loads(str(field))

class SQLiteDataManager(DataManagerBase):
    _config_schema = \
    """
    CREATE TABLE we_config (
        name    TEXT NOT NULL,
        value    BLOB,
        PRIMARY KEY (name)
    )
    """
    
    _state_schema = \
    """
    CREATE TABLE we_state (
        name    TEXT NOT NULL,
        value    BLOB,
        PRIMARY KEY (name)
    )
    """
    
    _stats_schema = \
    """ 
    CREATE TABLE we_stats (
        we_iter            UNSIGNED INTEGER NOT NULL,
        n_particles        UNSIGNED INTEGER NOT NULL,
        norm               REAL NOT NULL,
        cputime            REAL,
        walltime           REAL,
        PRIMARY KEY (we_iter)
    )
    """
    
    _data_schema = \
    """
    CREATE TABLE we_data (
        we_iter     UNSIGNED INTEGER NOT NULL,
        name        TEXT NOT NULL,
        value       BLOB NULL,
        PRIMARY KEY (we_iter, name)
    )
    """
        
    _segment_schema = \
    """
    CREATE TABLE segments (
        we_iter             UNSIGNED INTEGER NOT NULL,
        seg_id              UNSIGNED INTEGER NOT NULL,
        status              UNSIGNED INTEGER NOT NULL,
        p_parent_id         UNSIGNED INTEGER NULL,
        data_ref            TEXT NULL,
    
        weight              REAL NOT NULL,
    
        cputime             REAL NULL,
        walltime            REAL NULL,
        
        final_pcoord        BLOB NULL,
        supplementary_data  BLOB NULL,
        
        PRIMARY KEY (we_iter, seg_id)
    )
    """
    
    _seg_lineage_schema = \
    """
    CREATE TABLE segment_lineage (
        we_iter   UNSIGNED INTEGER NOT NULL,
        seg_id    UNSIGNED INTEGER NOT NULL,
        parent_id UNSIGNED INTEGER NOT NULL,
        UNIQUE    (we_iter, seg_id, parent_id)
    )
    """
        
    _tables = (_config_schema, _state_schema, _stats_schema, _data_schema,
               _segment_schema, _seg_lineage_schema)

    def __init__(self, source):
        source = os.path.normpath(os.path.abspath(source))
        dname = os.path.dirname(source)
        fname = os.path.basename(source) 
        DataManagerBase.__init__(self, source)
        
    def _get_conn(self):
        try:
            return self._conn
        except AttributeError:
            self._conn = sqlitedb.connect(self.source, LOCK_TIMEOUT,
                                          isolation_level = None)
            self._conn.row_factory = sqlitedb.Row
            return self._conn
            
    conn = property(_get_conn)
            
    def initialize(self):
        try:
            mtime = os.path.getmtime(self.source)
        except OSError, e:
            if e.errno == errno.ENOENT:
                pass
        else:
            raise OSError(errno.EEXIST, os.strerror(errno.EEXIST), self.source)
        
        conn = self.conn
        for table_schema in self._tables:
            conn.execute(table_schema)

    def store_state_object(self, name, obj):
        conn = self.conn
        curs = conn.cursor()
        curs.execute('''INSERT OR REPLACE INTO we_state VALUES (?, ?)''',
                     (name, pickle_to_blob(obj)))
        conn.commit()
        curs.close()
        
    def load_state_object(self, name):
        curs = self.conn.cursor()
        curs.execute('''SELECT value FROM we_state WHERE name=?''', (name,))
        return unpickle_from_blob(curs.fetchone()[0])
                        
    def save_config(self):
        curs = self.conn.cursor()
        curs.executemany('''INSERT OR REPLACE INTO we_config VALUES (?, ?)''',
                         ((k, pickle_to_blob(v)) for k, v in self.config.iteritems()))
        self.conn.commit()
        curs.close()
        
    def restore_config(self):
        curs = self.conn.cursor()
        curs.execute('''SELECT name, value FROM we_config''')
        self.config = dict((row[0], unpickle_from_blob(row[1]))
                           for row in curs)        
                        
    def is_propagation_complete(self):
        conn = self.conn
        curs = conn.cursor()
        curs.execute('''SELECT COUNT(*) FROM segments 
                        WHERE status != ?''', 
                        (Segment.SEG_STATUS_COMPLETE,) )
        n_incomplete = curs.fetchone()[0]
        curs.close()
        return bool(n_incomplete == 0)
    
    def _segment_from_row(self, row):
        seg =  Segment(we_iter = row['we_iter'],
                       seg_id = row['seg_id'],
                       status = row['status'],
                       weight = row['weight'],
                       p_parent_id = row['p_parent_id'],
                       data_ref = row['data_ref'],
                       walltime = row['walltime'],
                       cputime = row['cputime'],
                       supplementary_data = unpickle_from_blob(row['supplementary_data']),
                       final_pcoord = unpickle_from_blob(row['final_pcoord']))
        return seg
            
    def get_segment(self, we_iter, seg_id, load_parent = False):
        curs = self.conn.cursor()
        curs.execute('''SELECT * FROM segments 
                        WHERE we_iter = ? AND seg_id = ?''',
                     (we_iter, seg_id,))
        row = curs.fetchone()
        if row is None: 
            return None
        segment = self._segment_from_row(row)
        
        if load_parent and segment.p_parent_id is not None:
            curs.execute('''SELECT * FROM segments 
                              WHERE we_iter = ? AND seg_id = ?''',
                         (we_iter-1, segment.p_parent_id,))
            segment.parent = self._segment_from_row(curs.fetchone())
        return segment
    
    _update_segment_valid_fields = frozenset(('status', 'weight', 
                                              'data_ref', 'final_pcoord',  
                                              'cputime', 'walltime',
                                              'supplementary_data'))
    _update_segment_stmt = '''UPDATE segments SET status=?, weight=?,
                                data_ref=?, cputime=?, walltime=?,
                                supplementary_data=?, final_pcoord=?
                              WHERE we_iter=? AND seg_id=?'''
    def update_segment(self, segment, fields = None):
        fields = fields or {}
        for (k,v) in fields.iteritems():
            if k not in self._update_segment_valid_fields:
                raise TypeError('invalid field %r specified' % k)
            else:
                setattr(segment, k, v)
                
        if segment.supplementary_data:
            sdata = pickle_to_blob(segment.supplementary_data)
        else:
            sdata = None
            
        if segment.final_pcoord is None:
            fpc = None
        else:
            fpc = pickle_to_blob(segment.final_pcoord)
        
        curs = self.conn.cursor()
        curs.execute('BEGIN IMMEDIATE TRANSACTION')
        curs.execute(self._update_segment_stmt,
                     (segment.status, segment.weight, segment.data_ref,
                      segment.cputime, segment.walltime, sdata, fpc,
                      segment.we_iter, segment.seg_id))
        curs.execute('COMMIT')
        self.conn.commit()

    def get_segments(self, we_iter):
        curs = self.conn.cursor()
        curs.execute('''SELECT * FROM segments
                        WHERE we_iter = ?''', (we_iter,))
        segments = [self._segment_from_row(row) for row in curs]
        return segments
    
    def get_total_cputime(self):
        import datetime
        curs = self.conn.cursor()
        curs.execute('''SELECT SUM(cputime) FROM segments''')
        return datetime.timedelta(seconds = curs.fetchone()[0] or 0)
        
    def create_segments(self, segments):
        curs = self.conn.cursor()
        curs.execute('BEGIN IMMEDIATE TRANSACTION')
        curs.executemany('''INSERT INTO segments
                            (seg_id, we_iter, status, p_parent_id, weight,
                             data_ref, walltime, cputime, final_pcoord,
                             supplementary_data)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                         ((seg.seg_id, seg.we_iter, seg.status, seg.p_parent_id,
                           seg.weight,
                           seg.data_ref, seg.walltime, seg.cputime,
                           seg.final_pcoord is not None
                             and pickle_to_blob(seg.final_pcoord)
                             or None, 
                           seg.supplementary_data is not None
                             and pickle_to_blob(seg.supplementary_data)
                             or None) for seg in segments))
        curs.execute('COMMIT')
        self.conn.commit()
                
    _scrub_fields = ('cputime', 'walltime', 'final_pcoord', 'supplementary_data')
    def scrub_crashed_segments(self, we_iter):
        curs = self.conn.cursor()
        curs.execute('BEGIN EXCLUSIVE TRANSACTION')
        sstmt = ', '.join('%s=NULL' % field for field in self._scrub_fields)
        curs.execute('''UPDATE segments SET %s,
                                            status=?
                        WHERE we_iter=? AND status != ?''' % sstmt,
                     (Segment.SEG_STATUS_PREPARED, 
                      we_iter, Segment.SEG_STATUS_COMPLETE))
        curs.commit('COMMIT')
        self.conn.commit()
                
    def get_next_segment(self, pretend = False):
        curs = self.conn.cursor()
        if not pretend:
            curs.execute('''BEGIN EXCLUSIVE TRANSACTION''')
        
        curs.execute('''SELECT * FROM segments
                        WHERE we_iter = ? AND status = ?
                        ORDER BY seg_id LIMIT 1''',
                     (self.we_sim.current_iteration, 
                      Segment.SEG_STATUS_PREPARED))
        try:
            segment = self._segment_from_row(curs.fetchone())
        except (IndexError,TypeError):
            return None
        
        if segment.p_parent_id is not None:
            curs.execute('''SELECT * FROM segments 
                            WHERE we_iter = ? AND seg_id = ?''',
                         (segment.we_iter-1,segment.p_parent_id,))
            segment.p_parent = self._segment_from_row(curs.fetchone())
        
        if not pretend:
            curs.execute('''UPDATE segments SET status = ? 
                            WHERE we_iter = ? AND seg_id = ?''',
                            (Segment.SEG_STATUS_RUNNING, 
                             segment.we_iter, segment.seg_id))
        curs.execute('COMMIT')
        self.conn.commit()
        return segment
        
    def record_lineage(self, segment_lineages):
        curs = self.conn.cursor()
        curs.executemany('''INSERT INTO segment_lineage 
                              (we_iter, seg_id, parent_id)
                               VALUES (?, ?, ?)''',
                         ((segment.we_iter, segment.seg_id, parent.seg_id)
                          for (segment, parent) in segment_lineages))
        self.conn.commit()


    def record_data_item(self, we_iter, name, obj):
        curs = self.conn.cursor()
        curs.execute('''INSERT INTO we_data 
                          (we_iter, name, value) VALUES (?, ?, ?)''',
                     (we_iter, name, pickle_to_blob(obj)))
        
    def load_data_item(self, we_iter, name):
        curs = self.conn.cursor()
        curs.execute('SELECT value FROM we_data WHERE we_iter=? AND name=?',
                     (we_iter, name))
        row = curs.fetchone()
        if row is None:
            return None
        else:
            return unpickle_from_blob(row[0])
        
    def load_data_sequence(self, name):
        curs = self.conn.cursor()
        curs.execute('''SELECT we_iter, value FROM we_data WHERE name=?
                        ORDER BY we_iter ASC''', (name,))
        return [(row[0], unpickle_from_blob(row[1])) for row in curs]
