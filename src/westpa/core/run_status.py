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


@dataclass
class RunStatusReadResult:
    path: str
    status: dict | None = None
    error: str | None = None
    missing: bool = False


def status_path_for_h5(we_h5filename):
    return os.path.abspath(we_h5filename) + '.progress.json'


def read_run_status(we_h5filename):
    path = status_path_for_h5(we_h5filename)
    try:
        with open(path, 'r', encoding='utf-8') as status_file:
            status = json.load(status_file)
    except FileNotFoundError:
        return RunStatusReadResult(path=path, missing=True)
    except Exception as exc:
        return RunStatusReadResult(path=path, error=f'Could not read live status file {path}: {exc}')

    if not isinstance(status, dict):
        return RunStatusReadResult(path=path, error=f'Live status file {path} does not contain a JSON object')

    schema_version = status.get('schema_version')
    if schema_version != SCHEMA_VERSION:
        return RunStatusReadResult(
            path=path,
            error=f'Live status file {path} has unsupported schema version {schema_version!r}',
        )

    return RunStatusReadResult(path=path, status=status)


