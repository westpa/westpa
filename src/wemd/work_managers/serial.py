from wemd.work_managers import WorkManager
from wemd.core.errors import PropagationIncompleteError

class SerialWorkManager(WorkManager):
    def propagate_segments(self, we_iter):        
        segments = self.sim_driver.q_incomplete_segments(we_iter).all()
        dbsession = self.sim_driver.dbsession
        for segment in segments:
            self.log.info('dispatching segment %s (weight %g)' 
                          % (segment.seg_id, segment.weight))
            dbsession.begin()
            try:
                self.backend_driver.propagate_segment(segment)
                dbsession.flush()
            except:
                self.log.debug('error merging propagation results', exc_info = True)   
                dbsession.rollback()
                raise
            else:
                dbsession.commit()
