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


class Test_W_Progress:
    def test_snapshot_completed_run(self, ref_50iter):
        snapshot = read_progress_snapshot(self.h5_filepath, requested_total_iterations=50, recent=5)

        assert snapshot.error is None
        assert snapshot.current_iteration == 51
        assert snapshot.latest_completed_iteration == 50
        assert snapshot.requested_total_iterations == 50
        assert snapshot.completed_segments == 9985
        assert snapshot.completed_walltime > 0
        assert snapshot.average_recent_walltime > 0
        assert snapshot.eta_seconds == 0

    def test_snapshot_initialized_run(self, ref_initialized):
        snapshot = read_progress_snapshot(self.h5_filepath, requested_total_iterations=2, recent=5)

        assert snapshot.error is None
        assert snapshot.current_iteration == 1
        assert snapshot.latest_completed_iteration == 0
        assert snapshot.current_status_counts.total > 0
        assert snapshot.current_status_counts.complete == 0
        assert snapshot.current_status_counts.prepared == snapshot.current_status_counts.total
        assert snapshot.completed_segments is None
        assert snapshot.eta_seconds is None

    def test_missing_file_is_error_snapshot(self, ref_50iter):
        snapshot = read_progress_snapshot('missing-west.h5', requested_total_iterations=50, recent=5)
        output = render_progress(snapshot, refresh_interval=1.0)

        assert snapshot.error is not None
        assert 'Error:' in output
        assert 'missing-west.h5' in output

