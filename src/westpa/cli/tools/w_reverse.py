import logging
import os
import shutil
import tempfile

import numpy as np
from tqdm.auto import tqdm, trange

from westpa.core.h5io import WESTIterationFile
from westpa.core.data_manager import seg_id_dtype, n_iter_dtype, weight_dtype
from westpa.core.propagators.loaders import restart_writer
from westpa.core.segment import Segment
from westpa.core.trajectory import mdtraj_supported_extensions
from westpa.core.yamlcfg import YAMLConfig
from westpa.tools import WESTTool, WESTDataReader

log = logging.getLogger('w_reverse')


class W_Reverse(WESTTool):
    prog = 'w_reverse'
    description = '''
w_reverse: a tool for taking a WE simulation facilitated
through WESTPA with successful recycling events and generating
a new directory with the recycled restart files, which can then
serve as the bstates for a subsequent WE simulation in the opposite
direction, i.e. starting from the successfully recycled events and
then going back to the original starting point; overall reversed.

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output directory (--output-bstates-dir,-obd, by default "bstates_reverse") contains the following files:

  Trajectory restart files
    One trajectory file is saved for each successful trajctory segment in the simulation. If there are more than ``--max-n-bstates``
    successful recycled trajectory segments (default: 10000), then ``--max-n-bstates`` number of segments will be randomly picked
    and saved. The output files will be named ``{iteration:06d}_{walker:06d}`` followed by its original file extension.

  Bstate file
    File named ``--output-bstates-file``, by default ``bstates.txt``, contains the bstate index, weight, and file name for
    each copied trajectory file in the output directory.

'''

    def __init__(self):
        super().__init__()
        self.config = YAMLConfig()
        self.data_reader = WESTDataReader()
        self.rng = None
        self.top_exts, self.traj_exts = mdtraj_supported_extensions()

    @property
    def h5(self):
        return self.data_reader.data_manager.we_h5file

    def add_args(self, parser):
        self.data_reader.add_args(parser)

        rgroup = parser.add_argument_group('w_reverse options')
        rgroup.add_argument(
            '-W',
            '--west-data',
            dest='we_h5filename',
            metavar='WEST_H5FILE',
            type=str,
            default='west.h5',
            help='''Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in %(default)s).''',
        )
        rgroup.add_argument(
            '-r',
            '--rcfile',
            metavar='RCFILE',
            dest='rcfile',
            type=str,
            default='west.cfg',
            help='use RCFILE as the WEST run-time configuration file (default: %(default)s)',
        )
        rgroup.add_argument(
            "--first-iter", "-fi", dest="first_iter", type=int, default=1, help="First iteration to consider (default: %(default)s)"
        )
        rgroup.add_argument(
            "--last-iter",
            "-li",
            dest="last_iter",
            type=int,
            default=None,
            help="Last iteration to consider (default: None, i.e. last completed iteration in WEST_H5FILE)",
        )
        rgroup.add_argument(
            "--max-n-bstates",
            '--max',
            '-m',
            dest="max_n_bstates",
            type=int,
            default=10000,
            help="Max number of bstates to copy over. Adjust this if you prefer to include"
            + "a larger subset of the first bstates found (default: max of %(default)s).",
        )
        rgroup.add_argument("--rst-file", '-rf', dest="rst_file", type=str, default=None, help="Path to the Restart File")
        rgroup.add_argument(
            "--output-bstates-dir",
            "-od",
            "-o",
            dest="output_bstates_dir",
            type=str,
            default="bstates_reverse",
            help="Output directory for the bstates and output_bstates_file",
        )
        rgroup.add_argument(
            "--output-bstates-file",
            '-of',
            dest="output_bstates_file",
            type=str,
            default="bstates.txt",
            help="Name of the output bstates file",
        )
        rgroup.add_argument(
            "--use-weights",
            '--weights',
            "-uw",
            action='store_true',
            dest="use_weights",
            help="Include the recycled event weight when making the bstate file",
        )
        rgroup.add_argument(
            '--seed',
            '-s',
            type=int,
            default=None,
            dest='seed',
            help='Psuedo-random number generator seed used to select output structures.',
        )

    def process_args(self, args):
        # Process config from parser, rcfile ('west.cfg')
        self.data_reader.process_args(args)
        self.config.update_from_file(args.rcfile)
        self.data_refs_dic = self.config.get(['west', 'data', 'data_refs'], {})

        ## Open the west.h5 file
        # self.h5 = self.data_reader.data_manager.we_h5file

        # HDF5 framework or file path wrangling
        self.h5_framework = True if 'iteration' in self.data_refs_dic else False
        dict_key = 'iteration' if self.h5_framework else 'segment'
        self.traj_seg_path = os.path.expandvars(
            self.data_refs_dic[dict_key].replace('segment.n_iter', 'n_iter').replace('segment.seg_id', 'seg_id')
        )
        self.traj_seg_path = (
            self.traj_seg_path.replace('$WEST_SIM_ROOT', '.') if 'WEST_SIM_ROOT' in os.environ else self.traj_seg_path
        )
        self.output_bstates_dir = args.output_bstates_dir
        self.output_bstates_file = args.output_bstates_file

        # File extension wrangling
        self.rst_file = args.rst_file if args.rst_file else None
        self.rst_extension = self.rst_file.lower().rsplit('.', maxsplit=1)[-1] if self.rst_file else None
        self.traj_exc_exts = []
        self.traj_or_top_exts = []
        for i in self.traj_exts:
            if i in self.top_exts:
                self.traj_or_top_exts.append(i)
            else:
                self.traj_exc_exts.append(i)

        # Dealing with other w_reverse parameters
        with self.data_reader:
            self.first_iter = args.first_iter
            self.last_iter = args.last_iter or self.h5.attrs['west_current_iteration'] - 1
            self.max_n_bstates = args.max_n_bstates
            self.use_weights = args.use_weights

        # pRNG stuff
        self.seed = args.seed
        log.info(f'Using seed: {self.seed}')

    def _find_recycled_segments(self):
        """
        Find and return an array containing all successfully recycled segments' (iter_id, seg_id, weight)
        tuples based on the the segment endpoint status in the main HDF5 file (`west.h5`).

        Returns
        -------
        succ : np.ndarray of shape (n_successful_segments, 3)
            An array of all successfully recycled segments. Each row is [iteration, walker, weight].

        """
        succ = []
        for iteration in trange(
            1, len(self.h5['summary'][self.h5['summary']['walltime'] != 0]) + 1, desc="finding successful trajectories", leave=False
        ):
            endpoint_type = self.h5[f'iterations/iter_{iteration:>08d}/seg_index']['endpoint_type', :]
            indices = np.flatnonzero(endpoint_type == Segment.SEG_ENDPOINT_RECYCLED)
            temp_list = [
                (
                    iteration,
                    index,
                    self.h5[f'iterations/iter_{iteration:>08d}/seg_index']['weight', index],
                )
                for index in indices
            ]
            succ += temp_list

        dtypes = np.dtype([('n_iter', n_iter_dtype), ('seg_id', seg_id_dtype), ('weight', weight_dtype)])

        return np.asarray(succ, dtype=dtypes)

    def copy_traj(self, search_folder, iteration, walker):
        """
        Find a restart file in the ``search_folder`` and copying it to the ``output_bstates_dir`` under the name
        ``{iteration:06d}_{walker:06d}`` and the same file extension. Priority given if ``--rst-file`` is provided,
        first based on exact name matching or if not found, a file with the same extension with the newest creation time.
        Else, it will attempt to find the newest created file with a Mdtraj supported format extension.

        Arguments
        ---------
        search_folder: Path or str
            Directory path to be searched for a trajectory file.

        iteration: int
            Iteration number of the trajectory segment.

        walker: int
            Segment id of the trajectory segment.

        Returns
        -------
        rst_dest_name: str
            Name of the copied trajectory in the ``output_bstates_dir``

        """
        files = [file for file in os.listdir(search_folder) if not file.startswith('.')]
        if self.rst_file in files:
            # Exact file match
            source_file = self.rst_file
            rst_dest_name = f"{iteration:06d}_{walker:06d}.{self.rst_extension}"
        elif self.rst_extension:
            # Attempt to find based on file extension provided through `--rst-file` / `self.rst_file`
            possible_hits = [file for file in files if file.endswith(self.rst_extension)]
            if len(possible_hits) > 1:
                possible_hits = sorted(possible_hits, key=lambda file: os.path.getctime(file))
                log.warning(
                    f'Found {possible_hits[-1]} as restart file for iteration {iteration} and walker {walker} based on file creation times. if this is incorrect, provide a file name using flag --rst-file'
                )
            rst_dest_name = f"{iteration:06d}_{walker:06d}.{self.rst_extension}"
            source_file = possible_hits[-1]
            shutil.copyfile(os.path.join(search_folder, possible_hits[-1]), os.path.join(self.output_bstates_dir, rst_dest_name))
        else:
            # Guess based on MDTraj-supported file extensions
            possible_hits = []
            for test_extension in self.traj_exc_exts:
                possible_hits += [file for file in files if file.endswith(test_extension)]
            if len(possible_hits) < 1:
                for test_extension in self.traj_or_top_exts:
                    possible_hits += [file for file in files if file.endswith(test_extension)]
            elif len(possible_hits) >= 1:
                possible_hits = sorted(possible_hits, key=lambda file: os.path.getctime(os.path.join(search_folder, file)))
                log.warning(
                    f'Found {possible_hits[-1]} as restart file for iteration {iteration} and walker {walker} based on file creation times. if this is incorrect, provide a file name using flag `--rst-file`.'
                )
            source_file = possible_hits[-1]
            self.rst_extension = possible_hits[-1].split('.')[-1] if self.rst_extension is None else self.rst_extension
            rst_dest_name = f"{iteration:06d}_{walker:06d}.{possible_hits[-1].split('.')[-1]}"

        shutil.copyfile(os.path.join(search_folder, source_file), os.path.join(self.output_bstates_dir, rst_dest_name))

    def go(self):
        """
        Main public method for running w_reverse. First runs ``_find_recycled_segments()`` to find the successfully
        recycled trajectories, then iterates over them to copy the trajectories to ``output_bstates_dir``. Lastly,
        iterate over the trajectories again to make the ``output_bstates_file``.

        """
        self.rng = np.random.default_rng(seed=self.seed) if self.rng is None else self.rng

        with self.data_reader:

            # make directory for bstates_reverse if it doesn't already exist
            os.makedirs(self.output_bstates_dir, exist_ok=True)

            # Find successfully recycled trajectories and pick trajectories to copy
            succ_pairs = self._find_recycled_segments()
            total_pairs = min(self.max_n_bstates, len(succ_pairs))
            succ_pairs_used = self.rng.choice(
                succ_pairs,
                size=total_pairs,
                p=succ_pairs['weight'] / np.sum(succ_pairs['weight']),
                replace=False,
            )

        # Loop though picked segments and copy
        total_weight = 0.0
        for idx, succ_pair in enumerate(tqdm(succ_pairs_used, total=total_pairs, desc="New bstates")):
            [iteration, walker, weight] = succ_pair
            total_weight += weight
            # check if using HDF5 framework
            if self.h5_framework:
                with tempfile.TemporaryDirectory() as tmpdirname:
                    # Extract the restart data from the .h5 file
                    segment = Segment(n_iter=iteration, seg_id=walker, weight=weight)
                    with WESTIterationFile(self.traj_seg_path.format(n_iter=iteration)) as h5file:
                        h5file.read_restart(segment)
                        restart_writer(tmpdirname, segment)
                        self.copy_traj(tmpdirname, iteration, walker)
            else:
                search_folder = self.traj_seg_path.format(n_iter=iteration, seg_id=walker)
                self.copy_traj(search_folder, iteration, walker)

        # Loop through picked segemnts again and create `bstates.txt` file
        # fill out the bstates.txt file with name and weight
        # but only use weights if requested, otherwise use equal weights
        # bstates.txt row format: bstate_n | weight | bstate_filename
        with open(os.path.join(self.output_bstates_dir, self.output_bstates_file), "w") as bstates_f:
            for idx, succ_pair in enumerate(tqdm(succ_pairs_used, desc='writing bstate file', leave=False)):
                [iteration, walker, weight] = succ_pair
                rst_dest_name = f'{iteration:06d}_{walker:06d}.{self.rst_extension}'

                # Do some quick file checking to make sure files are there...
                files = [file for file in os.listdir(self.output_bstates_dir) if file.startswith(f'{iteration:06d}_{walker:06d}')]
                if len(files) > 1 and rst_dest_name not in files:
                    files = sorted(files, key=lambda file: os.path.getctime(file))
                    rst_dest_name = files[-1]
                    log.warning(
                        f'Muliple files starting with `{iteration:06d}_{walker:06d}` are in the output directory. Using {rst_dest_name=}. This is usually caused by running multiple rounds of `w_reverse` with different parameters.'
                    )
                else:
                    log.error(
                        f'A restart file starting with {iteration:06d}_{walker:06d} should be present in the output directory but is not!!!'
                    )

                # Actual writing into the file
                if self.use_weights:
                    bstates_f.write(f"{idx} {(weight / total_weight):.3e} {rst_dest_name}\n")
                else:
                    bstates_f.write(f"{idx} {(1 / total_pairs):.3e} {rst_dest_name}\n")


def entry_point():
    W_Reverse().main()


if __name__ == "__main__":
    entry_point()
