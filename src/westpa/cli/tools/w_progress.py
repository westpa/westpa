import sys
import time

import westpa
from westpa.core.run_status import ACTIVE_RUN_STATES, RUN_STATE_COMPLETE, read_run_status
from westpa.tools import WESTTool
from westpa.tools.progress_status import (
    progress_snapshot_from_run_status,
    read_progress_snapshot,
    render_progress,
)


class WProgress(WESTTool):
    prog = 'w_progress'
    description = 'Show a live progress dashboard for a WESTPA simulation.'

    def __init__(self):
        super().__init__()
        self.we_h5filename = None
        self.refresh_interval = 1.0
        self.requested_total_iterations = None

    def add_args(self, parser):
        group = parser.add_argument_group('WEST input data options')
        group.add_argument(
            '-W',
            '--west-data',
            dest='we_h5filename',
            metavar='WEST_H5FILE',
            help='Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in west.cfg).',
        )
        group.add_argument(
            '--refresh',
            dest='refresh_interval',
            type=float,
            default=1.0,
            metavar='SECONDS',
            help='Refresh the dashboard every SECONDS seconds (default: 1.0).',
        )

