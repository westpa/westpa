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
        * write a simple test: this should be very straitforward (assert time of test h5 file is as expected)
        * update docs
    """

    prog = "w_timings"
    description = "A tool for aggregate simulation and wallclock time extraction."

    def __init__(self, tau=100):
        """
        Parameters
        ----------
        tau : int
            WESTPA dynamics propagation time in picoseconds. Default 100 = 100ps.

        """
        super().__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection(self.data_reader)
        self.tau = tau

    def go(self):
        with self.data_reader:
            we_h5file = self.data_reader.data_manager.we_h5file
            iter_start = self.iter_range.iter_start
            iter_stop = self.iter_range.iter_stop

            walltime = we_h5file['summary']['walltime'][iter_start - 1 : iter_stop].sum()
            aggtime = we_h5file['summary']['n_particles'][iter_start - 1 : iter_stop].sum()

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

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.tau = args.tau
        with self.data_reader:
            self.iter_range.process_args(args)


def entry_point():
    WTimings().main()


if __name__ == "__main__":
    entry_point()
