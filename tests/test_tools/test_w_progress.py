import argparse
from unittest import mock

from westpa.cli.tools.w_progress import WProgress, entry_point
from westpa.core.run_status import RUN_STATE_COMPLETE, RUN_STATE_RUNNING, RunStatusReadResult
from westpa.tools.progress_status import (
    ProgressSnapshot,
    SegmentStatusCounts,
    format_duration,
    progress_snapshot_from_run_status,
    read_progress_snapshot,
    render_progress,
)


def live_status(**overrides):
    status = {
        'west_h5file': 'west.h5',
        'updated_at': 0,
        'run_state': RUN_STATE_RUNNING,
        'phase': 'propagating',
        'current_iteration': 37,
        'latest_completed_iteration': 36,
        'requested_total_iterations': 50,
        'segment_total': 100,
        'segment_prepared': 18,
        'segment_failed': 0,
        'iteration_started_at': None,
        'recent_walltimes': [],
        'completed_walltime': 642,
        'completed_segments': 9985,
        'message': None,
    }
    status.update(overrides)
    return status


