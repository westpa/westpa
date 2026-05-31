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


def _json_default(value):
    if hasattr(value, 'item'):
        return value.item()
    raise TypeError(f'Object of type {value.__class__.__name__} is not JSON serializable')


def write_run_status(we_h5filename, status):
    path = status_path_for_h5(we_h5filename)
    directory = os.path.dirname(path) or '.'
    os.makedirs(directory, exist_ok=True)

    payload = dict(status)
    payload['schema_version'] = SCHEMA_VERSION
    payload['west_h5file'] = os.path.abspath(we_h5filename)
    payload.setdefault('updated_at', time.time())

    tmp_name = None
    try:
        with tempfile.NamedTemporaryFile(
            'w',
            encoding='utf-8',
            dir=directory,
            prefix=os.path.basename(path) + '.',
            suffix='.tmp',
            delete=False,
        ) as tmp_file:
            tmp_name = tmp_file.name
            json.dump(payload, tmp_file, default=_json_default, sort_keys=True)
            tmp_file.write('\n')
        os.replace(tmp_name, path)
    finally:
        if tmp_name and os.path.exists(tmp_name):
            try:
                os.unlink(tmp_name)
            except OSError:
                pass


