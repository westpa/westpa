from wemd.work_managers import WorkManager
from wemd.core.errors import PropagationIncompleteError

import logging
log = logging.getLogger(__name__)

class SerialWorkManager(WorkManager):
    def propagate_segments(self, we_iter):
        n_inc = self.sim_driver.num_incomplete_segments(we_iter)
        log.info('%d segments remaining in WE iteration %d'
                 % (n_inc, we_iter))
        
        if n_inc == 0:
            return
        
        segments = self.sim_driver.get_incomplete_segments(we_iter)
        dbsession = self.sim_driver.dbsession
        for segment in segments:
            log.info('dispatching segment %s (weight %g)' 
                     % (segment.seg_id, segment.weight))
            dbsession.begin()
            try:
                self.backend_driver.propagate_segment(segment)
                dbsession.flush()
            except:
                log.debug('error merging propagation results', exc_info = True)   
                dbsession.rollback()
                raise
            else:
                dbsession.commit()
