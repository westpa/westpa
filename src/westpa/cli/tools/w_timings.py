import numpy as np
from tqdm.auto import tqdm

from westpa.tools import (
    WESTTool,
    WESTDataReader,
    IterRangeSelection,
)


class WTimings(WESTTool):
    """
    Aggregate simulation and wallclock time extraction.

    TODO:
        * convert to use logging?
        * convert from tqdm to native westpa progress indicator?
        * write a simple test: this should be very straitforward (assert time of test h5 file is as expected)
        * update docs
    """

    prog = "w_timings"
    description = "A tool for aggregate simulation and wallclock time extraction."

    def __init__(self, tau=100, count_events=False):
        """
        Parameters
        ----------
        tau : int
            WESTPA dynamics propagation time in picoseconds. Default 100 = 100ps.
        count_events : bool
            Option to also output the number of successfull recycling events.
        """
        super().__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection(self.data_reader)
        self.iter_start = None
        self.iter_stop = None
        self.tau = tau
        self.count_events = count_events

    def get_event_count(self, we_h5file):
        """
        Check if the target state was reached, given the data in a WEST H5 file.

        Parameters
        ----------
        we_h5file : h5py.File

        Returns
        -------
        int
            Number of successful recycling events.
        """
        events = 0
        # Get the key to the final iteration.
        # Need to do -2 instead of -1 because there's an empty-ish final iteration written.
        for iteration_key in tqdm(list(we_h5file['iterations'].keys())[-2:0:-1]):
            endpoint_types = we_h5file[f'iterations/{iteration_key}/seg_index']['endpoint_type']
            if 3 in endpoint_types:
                # print(f"recycled segment found in file {h5_filename} at iteration {iteration_key}")
                # count the number of 3s
                events += np.count_nonzero(endpoint_types == 3)
        return events

    def go(self):
        with self.data_reader:
            we_h5file = self.data_reader.data_manager.we_h5file
            self.iter_start = self.iter_range.iter_start
            self.iter_stop = self.iter_range.iter_stop

            walltime = we_h5file['summary']['walltime'][self.iter_start - 1 : self.iter_stop].sum()
            aggtime = we_h5file['summary']['n_particles'][self.iter_start - 1 : self.iter_stop].sum()

            days = int(walltime) // 86400
            hours = (int(walltime) % 86400) // 3600
            minutes = (int(walltime) % 3600) // 60
            seconds = walltime % 60

            print("\n===== WALLCLOCK  =====")
            print(f"{'Total Wallclock Time:':30}{days:>2}d {hours:>2}h {minutes:>2}m {seconds:>6.2f}s")

            print("\n===== SIMULATION  =====")
            print(f"{'Tau:':30} {self.tau} ps")
            print(f"{'Total Segments:':30} {aggtime}")
            print(f"{'Simulation time:':30} {((aggtime * self.tau) / 1000):.2f} ns")

            if self.count_events:
                events = self.get_event_count(we_h5file)
                print("\n===== RECYCLING =====")
                print(f"{'Recycled walkers:':30} {events}")

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
        parser.add_argument(
            "--tau",
            "-t",
            dest="tau",
            type=int,
            default=100,
            help="WESTPA dynamics propagation time in picoseconds. Default 100 = 100ps.",
        )
        parser.add_argument(
            "--count-events",
            "-ce",
            dest="count_events",
            action="store_true",
            default=False,
            help="Include this flag to also output the number of successfull recycling events.",
        )

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.tau = args.tau
        self.count_events = args.count_events
        with self.data_reader:
            self.iter_range.process_args(args)


def entry_point():
    WTimings().main()


if __name__ == "__main__":
    entry_point()
