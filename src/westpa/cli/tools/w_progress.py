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

