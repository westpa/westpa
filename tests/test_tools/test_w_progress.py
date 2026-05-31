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


