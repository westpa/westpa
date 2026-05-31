import os
import time
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np

from westpa.core.h5io import WESTPAH5File
from westpa.core.run_status import RUN_STATE_COMPLETE
from westpa.core.segment import Segment


@dataclass
class SegmentStatusCounts:
    total: int = 0
    unset: int = 0
    prepared: int = 0
    complete: int = 0
    failed: int = 0
    other: int = 0


@dataclass
class ProgressSnapshot:
    we_h5filename: str
    updated_at: float
    h5_mtime: float | None = None
    current_iteration: int | None = None
    latest_completed_iteration: int | None = None
    requested_total_iterations: int | None = None
    current_status_counts: SegmentStatusCounts = field(default_factory=SegmentStatusCounts)
    recent_walltimes: list[float] = field(default_factory=list)
    average_recent_walltime: float | None = None
    eta_seconds: float | None = None
    completed_walltime: float | None = None
    completed_segments: int | None = None
    error: str | None = None
    run_state: str | None = None
    phase: str | None = None
    status_updated_at: float | None = None
    current_iteration_elapsed: float | None = None
    message: str | None = None




def _current_iteration(h5file):
    attrs = h5file.attrs
    if 'west_current_iteration' in attrs:
        return int(attrs['west_current_iteration'])
    if 'wemd_current_iteration' in attrs:
        return int(attrs['wemd_current_iteration'])
    return None


def _segment_status_counts(iter_group):
    try:
        statuses = iter_group['seg_index']['status'][:]
    except KeyError:
        return SegmentStatusCounts()

    total = len(statuses)
    unset = int(np.count_nonzero(statuses == Segment.SEG_STATUS_UNSET))
    prepared = int(np.count_nonzero(statuses == Segment.SEG_STATUS_PREPARED))
    complete = int(np.count_nonzero(statuses == Segment.SEG_STATUS_COMPLETE))
    failed = int(np.count_nonzero(statuses == Segment.SEG_STATUS_FAILED))
    other = total - unset - prepared - complete - failed

    return SegmentStatusCounts(
        total=total,
        unset=unset,
        prepared=prepared,
        complete=complete,
        failed=failed,
        other=other,
    )


