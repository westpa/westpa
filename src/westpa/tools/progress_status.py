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


def progress_snapshot_from_run_status(status, requested_total_iterations=None, now=None):
    now = time.time() if now is None else now
    west_h5file = status.get('west_h5file', 'unknown')
    status_updated_at = status.get('updated_at')
    h5_mtime = None
    if west_h5file != 'unknown':
        try:
            h5_mtime = os.path.getmtime(west_h5file)
        except OSError:
            h5_mtime = None

    segment_total = int(status.get('segment_total') or 0)
    segment_prepared = int(status.get('segment_prepared') or 0)
    segment_failed = int(status.get('segment_failed') or 0)
    recent_walltimes = []
    for value in status.get('recent_walltimes', []):
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if _finite_positive(value):
            recent_walltimes.append(value)
    average_recent_walltime = None
    if recent_walltimes:
        average_recent_walltime = float(sum(recent_walltimes) / len(recent_walltimes))

    requested_total = status.get('requested_total_iterations', requested_total_iterations)
    if requested_total is not None:
        requested_total = int(requested_total)

    latest_completed = status.get('latest_completed_iteration')
    if latest_completed is not None:
        latest_completed = int(latest_completed)

    current_iteration = status.get('current_iteration')
    if current_iteration is not None:
        current_iteration = int(current_iteration)

    iteration_started_at = status.get('iteration_started_at')
    current_iteration_elapsed = None
    if iteration_started_at is not None:
        current_iteration_elapsed = max(now - float(iteration_started_at), 0.0)

    eta_seconds = _eta_seconds(requested_total, latest_completed, average_recent_walltime)
    if status.get('run_state') == RUN_STATE_COMPLETE:
        eta_seconds = 0

    return ProgressSnapshot(
        we_h5filename=west_h5file,
        updated_at=now,
        h5_mtime=h5_mtime,
        current_iteration=current_iteration,
        latest_completed_iteration=latest_completed,
        requested_total_iterations=requested_total,
        current_status_counts=SegmentStatusCounts(
            total=segment_total,
            prepared=segment_prepared,
            failed=segment_failed,
        ),
        recent_walltimes=recent_walltimes,
        average_recent_walltime=average_recent_walltime,
        eta_seconds=eta_seconds,
        completed_walltime=status.get('completed_walltime'),
        completed_segments=status.get('completed_segments'),
        run_state=status.get('run_state'),
        phase=status.get('phase'),
        status_updated_at=float(status_updated_at) if status_updated_at is not None else None,
        current_iteration_elapsed=current_iteration_elapsed,
        message=status.get('message'),
    )


def format_duration(seconds):
    if seconds is None or not np.isfinite(seconds):
        return 'unknown'
    seconds = float(seconds)
    if seconds < 60:
        return f'{seconds:.1f}s'

    seconds = int(round(seconds))
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f'{minutes}m {seconds:02d}s'

    hours, minutes = divmod(minutes, 60)
    if hours < 24:
        return f'{hours}h {minutes:02d}m {seconds:02d}s'

    days, hours = divmod(hours, 24)
    return f'{days}d {hours:02d}h {minutes:02d}m'


def _format_timestamp(timestamp):
    if timestamp is None:
        return 'unknown'
    return datetime.fromtimestamp(timestamp).strftime('%H:%M:%S')


def _format_value(value):
    return 'unknown' if value is None else str(value)


def _format_status_value(value):
    return 'unknown' if value is None else str(value).replace('_', ' ').title()


def _line(label, value, width=28):
    return f'{label:<{width}}{value}'


def _progress_text(latest_completed_iteration, requested_total_iterations):
    if latest_completed_iteration is None:
        return 'unknown'
    if requested_total_iterations is None or requested_total_iterations <= 0:
        return f'{latest_completed_iteration} completed'
    percent = 100.0 * latest_completed_iteration / requested_total_iterations
    return f'{latest_completed_iteration} / {requested_total_iterations} iterations ({percent:.1f}%)'


