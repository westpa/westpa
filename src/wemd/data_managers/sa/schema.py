__docformat__ = 'restructuredtext en'

import sqlalchemy
from sqlalchemy import (Table, Column, Index, ForeignKey,
                        SmallInteger, Integer, Boolean, Float, 
                        Text, CLOB, BLOB, PickleType, String, 
                        DateTime, Interval, UnicodeText)
from sqlalchemy.orm import (mapper, relation)
from sqlalchemy.orm.collections import column_mapped_collection
from sqlalchemy.ext.associationproxy import association_proxy

metadata = sqlalchemy.MetaData()

weIterTable = Table('we_iter', metadata,
                   Column('we_iter', Integer, primary_key=True),
                   Column('n_particles', Integer, nullable=False),
                   Column('norm', Float(17), nullable=False),
                   Column('cputime', Interval, nullable=True),
                   Column('walltime', Interval, nullable=True),
                   )

weDataTable = Table('we_data', metadata,
                  Column('we_iter', Integer, ForeignKey('we_iter.we_iter'),
                         primary_key=True, nullable=False,),
                  Column('name', Text, primary_key=True, nullable=False),
                  Column('pvalue', PickleType, nullable=True),
                  Column('cvalue', CLOB, nullable=True),
                  Column('bvalue', BLOB, nullable=True))

segDataTable = Table('seg_data', metadata,
                     Column('seg_id', Integer, ForeignKey('segments.seg_id'),
                            primary_key=True, nullable=False),
                     Column('name', Text, primary_key=True, nullable=False),
                     Column('pvalue', PickleType, nullable=True),
                     Column('cvalue', CLOB, nullable=True),
                     Column('bvalue', BLOB, nullable=True))

segmentsTable = Table('segments', metadata,
                      Column('seg_id', Integer, primary_key=True,
                             autoincrement=True),
                      Column('we_iter', Integer, nullable=False),                      
                      Column('status', SmallInteger, nullable=False),
                      Column('p_parent_id', Integer, 
                             ForeignKey('segments.seg_id'),
                             nullable=True),
                      Column('data_ref', Text, nullable=True),
                      Column('weight', Float(17), nullable=False),
                      Column('final_pcoord', PickleType, nullable=True),
                      Column('cputime', Interval, nullable=True),
                      Column('walltime', Interval, nullable=True),
                      Column('startdate', DateTime, nullable=True),
                      Column('enddate', DateTime, nullable=True),
                      )

Index('ix_segments_status', segmentsTable.c.we_iter,
                            segmentsTable.c.status)
Index('ix_segments_p_parent_id', segmentsTable.c.p_parent_id)

segmentLineageTable = Table('segment_lineage', metadata,
                            Column('seg_id', Integer,
                                   ForeignKey('segments.seg_id'), 
                                   primary_key=True, nullable=False),
                            Column('parent_id', Integer,
                                   ForeignKey('segments.seg_id'), 
                                   primary_key=True, nullable=False))

from wemd.core.segments import Segment
from wemd.core.we_sim import WESimIter
from data_items import DBDataItem

class SegmentDataItem(DBDataItem): pass
mapper(SegmentDataItem, segDataTable)

class WEIterDataItem(DBDataItem): pass
mapper(WEIterDataItem, weDataTable)

seg_pparent_rel = relation(Segment, segmentsTable, uselist=False)
seg_parents_rel = relation(Segment, segmentLineageTable,
                           collection_class = set,
                           primaryjoin=segmentsTable.c.seg_id==segmentLineageTable.c.seg_id,
                           secondaryjoin=segmentLineageTable.c.parent_id==segmentsTable.c.seg_id,
                           )
# Note that there is no delete-cascade if a parent gets deleted
# This could leave segmentLineageTable in an inconsistent state
# However, this shouldn't be trouble in practice because deleting a segment
# should only happen to delete a corrupted/incomplete segment in the current
# iteration, and parents are always in the previous iteration 

Segment.data = association_proxy('seg_data', 'value', creator=SegmentDataItem.named_create)
WESimIter.data = association_proxy('iter_data', 'value', creator=WEIterDataItem.named_create)

mapper(Segment, segmentsTable,
       properties = {'p_parent': seg_pparent_rel,
                     'parents': seg_parents_rel,
                     'seg_data': relation(SegmentDataItem, 
                                          collection_class=column_mapped_collection(segDataTable.c.name),
                                          cascade='all, delete, delete-orphan'), })

mapper(WESimIter, weIterTable,
       properties = {'iter_data': relation(WEIterDataItem,
                                           collection_class=column_mapped_collection(weDataTable.c.name),
                                           cascade='all, delete, delete-orphan'),})
