from wemd.work_managers import WorkManager

import logging
log = logging.getLogger(__name__)

class SerialWorkManager(WorkManager):
    def propagate_segments(self, segments):
        for segment in segments:
            log.info('dispatching segment %s' % segment.seg_id)
            self.backend_driver.propagate_segment(segment)
