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

    def process_args(self, args):
        if args.refresh_interval <= 0:
            self.parser.error('--refresh must be greater than 0')

        data_manager = westpa.rc.get_data_manager()
        if args.we_h5filename:
            data_manager.we_h5filename = args.we_h5filename

        self.we_h5filename = data_manager.we_h5filename
        self.refresh_interval = args.refresh_interval

        requested_total = westpa.rc.config.get(['west', 'propagation', 'max_total_iterations'], None)
        self.requested_total_iterations = int(requested_total) if requested_total is not None else None

    @property
    def should_clear(self):
        return sys.stdout.isatty()

    def snapshot(self):
        return read_progress_snapshot(
            self.we_h5filename,
            requested_total_iterations=self.requested_total_iterations,
            recent=5,
        )

    def sidecar_snapshot(self):
        status_result = read_run_status(self.we_h5filename)
        if status_result.status is None:
            return None, status_result
        return (
            progress_snapshot_from_run_status(
                status_result.status,
                requested_total_iterations=self.requested_total_iterations,
            ),
            status_result,
        )

    def render_once(self, include_hint=True):
        sidecar_snapshot, sidecar_status = self.sidecar_snapshot()
        if sidecar_snapshot is not None and sidecar_snapshot.run_state in ACTIVE_RUN_STATES:
            return render_progress(sidecar_snapshot, refresh_interval=self.refresh_interval, include_hint=include_hint)

        snapshot = self.snapshot()
        if snapshot.error is None:
            return render_progress(snapshot, refresh_interval=self.refresh_interval, include_hint=include_hint)

        if sidecar_snapshot is not None and sidecar_snapshot.run_state == RUN_STATE_COMPLETE:
            return render_progress(
                sidecar_snapshot,
                refresh_interval=self.refresh_interval,
                include_hint=include_hint,
                status_message=f'Could not read west.h5 ({snapshot.error}); showing the last live status snapshot instead.',
            )

        if sidecar_status.error:
            return render_progress(
                snapshot,
                refresh_interval=self.refresh_interval,
                include_hint=include_hint,
                status_message=sidecar_status.error,
            )

        return render_progress(snapshot, refresh_interval=self.refresh_interval, include_hint=include_hint)

    def go(self):
        try:
            while True:
                if self.should_clear:
                    sys.stdout.write('\033[H\033[J')
                sys.stdout.write(self.render_once())
                sys.stdout.flush()
                time.sleep(self.refresh_interval)
        except KeyboardInterrupt:
            if self.should_clear:
                sys.stdout.write('\n')
                sys.stdout.flush()


def entry_point():
    WProgress().main()


if __name__ == '__main__':
    entry_point()
