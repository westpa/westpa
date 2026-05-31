"""Live run status sidecar support for ``w_progress``.

The running ``w_run`` process owns ``west.h5`` while it is writing simulation
data, so ``w_progress`` uses this lightweight JSON file for active-run
monitoring and falls back to HDF5 for completed-run details.
"""

import json
import os
import tempfile
import time
from dataclasses import dataclass


SCHEMA_VERSION = 1
RUN_STATE_RUNNING = 'running'
RUN_STATE_COMPLETE = 'complete'
RUN_STATE_INTERRUPTED = 'interrupted'
RUN_STATE_ERROR = 'error'
ACTIVE_RUN_STATES = {RUN_STATE_RUNNING, RUN_STATE_INTERRUPTED, RUN_STATE_ERROR}


