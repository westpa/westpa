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


def _finite_positive(value):
    return value is not None and np.isfinite(value) and value > 0


def _h5_iter_group(h5file, n_iter):
    try:
        return h5file.get_iter_group(n_iter)
    except KeyError:
        return h5file[f'/iter_{int(n_iter):0{h5file.iter_prec}d}']


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


def _latest_completed_iteration(current_iteration, status_counts):
    if current_iteration is None:
        return None
    if status_counts.total and status_counts.complete == status_counts.total:
        return current_iteration
    return max(current_iteration - 1, 0)


def _summary_rows(h5file, latest_completed_iteration):
    if latest_completed_iteration is None or latest_completed_iteration <= 0:
        return None
    try:
        summary = h5file['summary']
    except KeyError:
        return None
    stop = min(latest_completed_iteration, len(summary))
    if stop <= 0:
        return None
    return summary[:stop]


def _walltimes(rows, recent):
    if rows is None or 'walltime' not in rows.dtype.names:
        return []
    values = [float(value) for value in rows['walltime'] if _finite_positive(float(value))]
    return values[-recent:]


def _completed_walltime(rows):
    if rows is None or 'walltime' not in rows.dtype.names:
        return None
    values = [float(value) for value in rows['walltime'] if np.isfinite(float(value))]
    return float(sum(values))


def _completed_segments(rows):
    if rows is None or 'n_particles' not in rows.dtype.names:
        return None
    return int(rows['n_particles'].sum())


def _eta_seconds(requested_total_iterations, latest_completed_iteration, average_recent_walltime):
    if requested_total_iterations is None or latest_completed_iteration is None:
        return None
    if not _finite_positive(average_recent_walltime):
        return None
    remaining_iterations = max(requested_total_iterations - latest_completed_iteration, 0)
    return remaining_iterations * average_recent_walltime


def read_progress_snapshot(we_h5filename, requested_total_iterations=None, recent=5, now=None):
    now = time.time() if now is None else now
    snapshot = ProgressSnapshot(
        we_h5filename=we_h5filename,
        requested_total_iterations=requested_total_iterations,
        updated_at=now,
    )

    try:
        snapshot.h5_mtime = os.path.getmtime(we_h5filename)
        with WESTPAH5File(we_h5filename, 'r') as h5file:
            current_iteration = _current_iteration(h5file)
            snapshot.current_iteration = current_iteration

            status_counts = SegmentStatusCounts()
            if current_iteration is not None and current_iteration > 0:
                try:
                    status_counts = _segment_status_counts(_h5_iter_group(h5file, current_iteration))
                except KeyError:
                    status_counts = SegmentStatusCounts()
            snapshot.current_status_counts = status_counts

            latest_completed = _latest_completed_iteration(current_iteration, status_counts)
            snapshot.latest_completed_iteration = latest_completed

            rows = _summary_rows(h5file, latest_completed)
            snapshot.recent_walltimes = _walltimes(rows, recent)
            if snapshot.recent_walltimes:
                snapshot.average_recent_walltime = float(sum(snapshot.recent_walltimes) / len(snapshot.recent_walltimes))
            snapshot.completed_walltime = _completed_walltime(rows)
            snapshot.completed_segments = _completed_segments(rows)
            snapshot.eta_seconds = _eta_seconds(
                requested_total_iterations,
                latest_completed,
                snapshot.average_recent_walltime,
            )
    except Exception as exc:
        snapshot.error = str(exc)

    return snapshot


