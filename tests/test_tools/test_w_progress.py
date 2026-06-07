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

    def test_render_progress(self):
        snapshot = ProgressSnapshot(
            we_h5filename='west.h5',
            updated_at=0,
            h5_mtime=0,
            current_iteration=37,
            latest_completed_iteration=36,
            requested_total_iterations=50,
            current_status_counts=SegmentStatusCounts(total=100, prepared=18, failed=2),
            recent_walltimes=[12.4, 11.8, 13.1],
            average_recent_walltime=12.433333333333334,
            eta_seconds=174.06666666666666,
            completed_walltime=642,
            completed_segments=9985,
        )

        output = render_progress(snapshot, refresh_interval=1.0)

        assert 'WESTPA Progress' in output
        assert 'Current iteration:' in output
        assert '37' in output
        assert '36 / 50 iterations (72.0%)' in output
        assert '18 / 100' in output
        assert '2 / 100' in output
        assert 'Recent walltimes:' in output
        assert '12.4s, 11.8s, 13.1s' in output
        assert 'ETA:' in output
        assert '2m 54s' in output

    def test_format_duration(self):
        assert format_duration(None) == 'unknown'
        assert format_duration(12.4) == '12.4s'
        assert format_duration(174.0) == '2m 54s'
        assert format_duration(3723.0) == '1h 02m 03s'

    def test_render_status_message(self):
        snapshot = ProgressSnapshot(
            we_h5filename='west.h5',
            updated_at=0,
            h5_mtime=0,
            current_iteration=1,
            latest_completed_iteration=0,
        )

        output = render_progress(snapshot, status_message='Using live status data.')

        assert 'Status' in output
        assert 'Using live status data.' in output

    def test_render_live_status_snapshot(self):
        snapshot = progress_snapshot_from_run_status(
            live_status(iteration_started_at=0, recent_walltimes=[12.4, 11.8, 13.1]),
            now=20,
        )

        assert snapshot.updated_at == 20
        assert snapshot.status_updated_at == 0
        output = render_progress(snapshot, refresh_interval=1.0)

        assert 'Run state:                  Running' in output
        assert 'Phase:                      Propagating' in output
        assert 'Live status updated:' in output
        assert 'Current iteration:          37' in output
        assert 'Prepared:                   18 / 100' in output
        assert 'Failed:                     0 / 100' in output
        assert 'Current iter elapsed:       20.0s' in output
        assert 'Completed walltime:         10m 42s' in output
        assert 'Completed segments:         9985' in output

    def test_parser(self, ref_50iter):
        tool = WProgress()
        tool.make_parser_and_process(args=['-W', self.h5_filepath, '--refresh', '2.5'])

        assert tool.we_h5filename == self.h5_filepath
        assert tool.refresh_interval == 2.5
        assert tool.requested_total_iterations == 50

    def test_render_once_running_sidecar_does_not_open_hdf5(self, ref_50iter):
        tool = WProgress()
        tool.we_h5filename = self.h5_filepath
        live_snapshot = progress_snapshot_from_run_status(
            live_status(west_h5file=self.h5_filepath)
        )

        with mock.patch.object(tool, 'sidecar_snapshot', return_value=(live_snapshot, RunStatusReadResult(path='status'))):
            with mock.patch.object(tool, 'snapshot', side_effect=AssertionError('west.h5 should not be opened')):
                output = tool.render_once(include_hint=False)

        assert 'Run state:                  Running' in output
        assert 'Current iteration:          37' in output

    def test_render_once_complete_sidecar_prefers_readable_hdf5(self, ref_50iter):
        tool = WProgress()
        tool.we_h5filename = self.h5_filepath
        sidecar_snapshot = progress_snapshot_from_run_status(
            live_status(west_h5file=self.h5_filepath, run_state=RUN_STATE_COMPLETE, phase='complete', segment_prepared=0)
        )
        hdf5_snapshot = ProgressSnapshot(
            we_h5filename=self.h5_filepath,
            updated_at=1,
            current_iteration=51,
            latest_completed_iteration=50,
            requested_total_iterations=50,
        )

        with mock.patch.object(tool, 'sidecar_snapshot', return_value=(sidecar_snapshot, RunStatusReadResult(path='status'))):
            with mock.patch.object(tool, 'snapshot', return_value=hdf5_snapshot):
                output = tool.render_once(include_hint=False)

        assert 'Current iteration:          51' in output
        assert 'Run state:' not in output

    def test_render_once_complete_sidecar_fallback_when_hdf5_unreadable(self, ref_50iter):
        tool = WProgress()
        tool.we_h5filename = self.h5_filepath
        sidecar_snapshot = progress_snapshot_from_run_status(
            live_status(
                west_h5file=self.h5_filepath,
                run_state=RUN_STATE_COMPLETE,
                phase='complete',
                current_iteration=51,
                latest_completed_iteration=50,
                segment_prepared=0,
                recent_walltimes=[12.0],
            )
        )
        unreadable_hdf5_snapshot = ProgressSnapshot(
            we_h5filename=self.h5_filepath,
            updated_at=1,
            error='Unable to open west.h5',
        )

        with mock.patch.object(tool, 'sidecar_snapshot', return_value=(sidecar_snapshot, RunStatusReadResult(path='status'))):
            with mock.patch.object(tool, 'snapshot', return_value=unreadable_hdf5_snapshot):
                output = tool.render_once(include_hint=False)

        assert 'Could not read west.h5' in output
        assert 'Run state:                  Complete' in output
        assert 'Error:' not in output

