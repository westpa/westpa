from datetime import timedelta

import numpy as np

from westpa.tools import (
    WESTTool,
    WESTDataReader,
    IterRangeSelection,
)


# NumPy fixed time units (excludes nonlinear units 'Y and 'M')
TIME_UNITS = ('W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns', 'ps', 'fs', 'as')


def _unit(delta):
    # Return the unit of a NumPy timedelta.
    words = str(delta.dtype).split('[')
    if len(words) == 1:
        return None
    return words[1][:-1]


def _str(delta):
    # Return a string representation of a NumPy timedelta.
    # Example: timedelta64(5380000, 'ps') -> '5.38 us'
    if _unit(delta) is None:
        return str(delta)
    for unit in TIME_UNITS:
        unit_delta = np.timedelta64(1, unit)
        try:
            if delta >= unit_delta:
                break
        except OverflowError:
            continue
    return f'{delta / unit_delta} {unit}'


class WTimings(WESTTool):
    prog = 'w_timings'
    description = 'Print timing information for a WESTPA simulation.'

    def __init__(self):
        super().__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection(self.data_reader)
        self.tau = None

    def go(self):
        start = self.iter_range.iter_start - 1
        stop = self.iter_range.iter_stop - 1
        with self.data_reader:
            iter_summaries = self.data_reader.we_h5file['summary'][start:stop]

        walltime = iter_summaries['walltime'].sum()
        cputime = iter_summaries['cputime'].sum()
        n_particles = iter_summaries['n_particles'].sum()
        n_iters = len(iter_summaries)

        print(f'Iterations: {n_iters}')
        print(f'Total segments: {n_particles}')
        print(f'Wall-clock time: {timedelta(seconds=walltime)}')
        if not np.isclose(cputime, 0):  # only print CPU time if it was recorded
            print(f'Total CPU time: {timedelta(seconds=cputime)}')
        if self.tau is not None:
            print(f'Simulated physical time ("molecular time"): {_str(n_iters * self.tau)}')
            print(f'Aggregate simulation time: {_str(n_particles * self.tau)}')

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)

        unit_list = ', '.join(map(repr, TIME_UNITS[-1:0:-1])) + f', or {TIME_UNITS[0]!r}'
        parser.add_argument(
            '-t',
            '--tau',
            help=(
                'WE resampling interval (format: <value>_<unit>, where <value> '
                f'is a positive integer and <unit> is {unit_list}).'
            ),
        )

    def process_args(self, args):
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)

        if args.tau is not None:
            value, unit = args.tau.split('_')
            if unit not in TIME_UNITS:
                raise ValueError(f'{unit!r} is not a recognized time unit')
            self.tau = np.timedelta64(int(value), unit)


def entry_point():
    WTimings().main()


if __name__ == "__main__":
    entry_point()
