'''Live progress dashboard for a running WESTPA simulation.'''

import os
import sys
import time
from datetime import datetime, timedelta

import numpy as np

import westpa
from westpa.core.h5io import WESTPAH5File
from westpa.core.segment import Segment
from westpa.tools import WESTTool, WESTDataReader


def read_progress(we_h5filename, max_total_iterations=None, n_recent=5):
    '''Return a dict of progress info from ``we_h5filename``.

    ETA is iterations left times the mean walltime of the last ``n_recent``
    iterations. None if ``max_total_iterations`` isn't set or nothing has
    finished yet.'''

    # w_run keeps west.h5 locked while it runs. Open without locking and read
    # whatever was last flushed.
    with WESTPAH5File(we_h5filename, 'r', locking=False) as h5file:
        n_iter = int(h5file.attrs['west_current_iteration'])
        n_completed = max(n_iter - 1, 0)
        summary = h5file['summary'][:n_completed]
        try:
            seg_status = h5file.get_iter_group(n_iter)['seg_index']['status']
        except KeyError:
            seg_status = np.empty((0,), dtype=np.uint8)

    # Ignore NaN and zero walltimes (truncated runs can leave these)
    walltimes = summary['walltime']
    recent_walltimes = walltimes[np.isfinite(walltimes) & (walltimes > 0)][-n_recent:]
    avg_walltime = float(recent_walltimes.mean()) if len(recent_walltimes) else None

    eta = None
    if max_total_iterations and avg_walltime is not None:
        eta = max(max_total_iterations - n_completed, 0) * avg_walltime

    return {
        'mtime': os.path.getmtime(we_h5filename),
        'n_iter': n_iter,
        'n_completed': n_completed,
        'max_total_iterations': max_total_iterations,
        'n_segs': len(seg_status),
        'n_complete': int(np.count_nonzero(seg_status == Segment.SEG_STATUS_COMPLETE)),
        'n_failed': int(np.count_nonzero(seg_status == Segment.SEG_STATUS_FAILED)),
        'n_particles': int(summary['n_particles'].sum()),
        'walltime': float(np.nansum(walltimes)),
        'recent_walltimes': recent_walltimes.tolist(),
        'avg_walltime': avg_walltime,
        'eta': eta,
    }


def _duration(seconds):
    '''Format seconds as H:MM:SS, or 'unknown' for None.'''
    if seconds is None:
        return 'unknown'
    return str(timedelta(seconds=round(seconds)))


def format_progress(progress):
    '''Build the dashboard text from a read_progress() dict.'''
    n_completed = progress['n_completed']
    max_total_iterations = progress['max_total_iterations']
    if max_total_iterations:
        percent = 100 * n_completed / max_total_iterations
        completed = f'{n_completed} / {max_total_iterations} iterations ({percent:.1f}%)'
    else:
        completed = f'{n_completed} completed'

    rows = [
        ('Last written:', datetime.fromtimestamp(progress['mtime']).strftime('%Y-%m-%d %H:%M:%S')),
        ('Current iteration:', progress['n_iter']),
        ('Progress:', completed),
        ('Segments complete:', f"{progress['n_complete']} / {progress['n_segs']}"),
        ('Segments failed:', progress['n_failed']),
        ('Recent walltimes:', ', '.join(map(_duration, progress['recent_walltimes'])) or 'unknown'),
        ('Avg iteration time:', _duration(progress['avg_walltime'])),
        ('Completed walltime:', _duration(progress['walltime'])),
        ('Completed segments:', progress['n_particles']),
        ('ETA:', _duration(progress['eta'])),
    ]

    width = 26
    return ''.join(f'{label.ljust(width)} {value}\n' for label, value in rows)


class WProgress(WESTTool):
    '''Redraw a progress dashboard until Ctrl-C.'''

    prog = 'w_progress'
    description = 'Show a live progress dashboard for a WESTPA simulation.'

    def __init__(self):
        super().__init__()
        self.data_reader = WESTDataReader()
        self.max_total_iterations = None
        self.refresh = None

    def add_args(self, parser):
        '''Add the -W and --refresh options.'''
        self.data_reader.add_args(parser)
        parser.add_argument(
            '--refresh',
            type=float,
            default=1.0,
            metavar='SECONDS',
            help='Redraw the dashboard every SECONDS seconds (default: %(default)s).',
        )

    def process_args(self, args):
        '''Get the HDF5 file and read max_total_iterations from west.cfg.'''
        self.data_reader.process_args(args)
        self.max_total_iterations = westpa.rc.config.get(['west', 'propagation', 'max_total_iterations'])

        if args.refresh <= 0:
            self.parser.error('argument --refresh: must be greater than 0')
        self.refresh = args.refresh

    def go(self):
        '''Redraw the dashboard every --refresh seconds.'''
        we_h5filename = self.data_reader.we_h5filename
        body = ''
        try:
            while True:
                try:
                    body = format_progress(read_progress(we_h5filename, self.max_total_iterations))
                    error = ''
                except Exception as e:
                    # A read can fail mid flush. Keep the last output and retry.
                    error = f'Could not read {we_h5filename}: {e}\n'

                if sys.stdout.isatty():
                    print('\033[H\033[J', end='')  # clear the terminal
                print(f'WESTPA progress for {we_h5filename} (updated {datetime.now():%H:%M:%S})\n')
                print(body + error, end='', flush=True)
                time.sleep(self.refresh)
        except KeyboardInterrupt:
            print()


def entry_point():
    WProgress().main()


if __name__ == '__main__':
    entry_point()
