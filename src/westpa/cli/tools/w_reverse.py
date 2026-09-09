import logging
from tqdm.auto import tqdm
import os
import shutil
import tempfile
from westpa.core.h5io import WESTIterationFile, WESTPAH5File
from westpa.core.propagators.loaders import restart_writer
from westpa.core.segment import Segment
from westpa.core._rc import WESTRC
from westpa.tools import WESTTool, WESTDataReader
import numpy as np
from westpa.core.trajectory import mdtraj_supported_extensions

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

  Trajectory files
    One trajectory file is saved for each successful trajctory in the simulation. If there are more than --max-n-bstates successful trajectories, by default 10000, then only --max-n-bstates trajectories will be saved in a random order. Trajectory files will be named {iteration:06d}_{walker:06d} with its original extension.

  Bstate file
    File named --output-bstates-file, by default bstates.txt, contains the index, weight, and file name for each trajectory file in the output directory

'''

    def __init__(self):
        super().__init__()
        self.westrc = WESTRC()
        self.data_reader = WESTDataReader()
        self.top_exts, self.traj_exts = mdtraj_supported_extensions()

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        rgroup = parser.add_argument_group('reverse options')
        rgroup.add_argument(
            "-W",
            "-w",
            "--west",
            "--west-data",
            "-h5",
            "--h5file",
            dest="we_h5filename",
            type=str,
            default="west.h5",
            help="Path to west.h5 file",
        )
        rgroup.add_argument(
            "--first-iter", "-fi", dest="first_iter", type=int, default=1, help="First iteration to consider (default: 1)"
        )
        rgroup.add_argument(
            "--last-iter",
            "-li",
            dest="last_iter",
            type=int,
            default=None,
            help="Last iteration to consider (default: last recorded iteration in west.h5)",
        )
        rgroup.add_argument("--config-file", dest="config_file", type=str, default="west.cfg", help="Path to the config file")
        rgroup.add_argument(
            "--max-n-bstates",
            dest="max_n_bstates",
            type=int,
            default=10000,
            help="Max number of bstates to copy over. Adjust this if you prefer "
            + "a subset of the first bstates found. Default max of 10000.",
        )
        rgroup.add_argument("--rst-file", dest="rst_file", type=str, default=None, help="Path to the Restart File")
        rgroup.add_argument(
            "--output-bstates-dir",
            "-obd",
            dest="output_bstates_dir",
            type=str,
            default="bstates_reverse",
            help="Output directory for the bstates and output_bstates_file",
        )
        rgroup.add_argument(
            "--output-bstates-file",
            "-obf",
            dest="output_bstates_file",
            type=str,
            default="bstates.txt",
            help="Name of the output bstates file",
        )
        rgroup.add_argument(
            "--use-weights",
            "-nw",
            type=bool,
            action='store_false',
            dest="use_weights",
            help="Include the recycled event weight when making the bstates.txt file",
        )
        rgroup.add_argument(
            '--seed', type=int, default=None, dest='seed', help='Seed for randomly choosing which states to include'
        )

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.config_required = True
        self.config_file = args.config_file
        self.westrc.read_config(self.config_file)
        self.config = self.westrc.config
        # Read the west.h5 file
        self.h5 = WESTPAH5File(args.we_h5filename, 'r')
        self.first_iter = args.first_iter
        self.last_iter = args.last_iter or self.h5.attrs['west_current_iteration'] - 1
        # Look at the data_refs from the config file
        self.data_refs_dic = self.config.get(['west', 'data', 'data_refs'], {})
        # Default to not using HDF5 framework
        self.h5_framework = True if 'iteration' in self.data_refs_dic else False
        dict_key = 'iteration' if self.h5_framework else 'segment'
        self.traj_seg_path = os.path.expandvars(
            self.data_refs_dic[dict_key].replace('segment.n_iter', 'n_iter').replace('segment.seg_id', 'seg_id')
        )
        # Proves that WEST_SIM_ROOT is not on os.environ, not replaced
        self.traj_seg_path = self.traj_seg_path.replace('$WEST_SIM_ROOT', '.')
        self.max_n_bstates = args.max_n_bstates
        self.rst_file = args.rst_file.lower() if args.rst_file else None
        self.rst_extension = self.rst_file.rsplit('.', maxsplit=1)[-1] if self.rst_file else None
        self.traj_exc_exts = []
        self.traj_or_top_exts = []
        for i in self.traj_exts:
            if i in self.top_exts:
                self.traj_or_top_exts.append(i)
            else:
                self.traj_exc_exts.append(i)
        self.output_bstates_dir = args.output_bstates_dir
        self.output_bstates_file = args.output_bstates_file
        self.use_weights = args.use_weights
        self.seed = args.seed
        log.info(f'Using seed: {self.seed}')

    def w_succ(self):
        """
        Find and return all successfully recycled (iter, seg) pairs.

        Returns
        -------
        succ : array of shape (n, 3) with [iteration, walker, weight] for each succ[i]
        """
        succ = []
        for iteration_index, iteration in tqdm(
            enumerate(self.h5['iterations'].keys()), total=len(self.h5['iterations'].keys()), desc="w_succ"
        ):
            endpoint_type = self.h5[f'iterations/{iteration}/seg_index']['endpoint_type', :]
            indices = np.flatnonzero(endpoint_type == Segment.SEG_ENDPOINT_RECYCLED)
            temp_array = [
                [
                    iteration_index if self.h5_framework else iteration_index + 1,
                    index,
                    self.h5[f'iterations/{iteration}/'seg_index']['weight', index],
                ]
                for index in indices
            ]
            succ += temp_array
        return np.asarray(succ)

    def copy_traj(self, search_folder, iteration, walker):
        """
        Find the correct trajectory file in the serch_folder and copy it to the output_bstates_dir with the name {iteration:06d}_{walker:06d} and keeping the same extension

        Arguments
        ---------
        search_folder: path or string describing a path to the directory that will be searched for a trajectory file

        iteration: integer describing the iteration number of the trajectory

        walker: integer describing the segment id of the trajectory

        Returns
        -------
        rst_dest_name: string describing the name of the copied trajectory in the output_bstates_dir

        """
        files = [file for file in os.listdir(search_folder) if not file.startswith('.')]
        if self.rst_file in files:
            rst_dest_name = f"{iteration:06d}_{walker:06d}.{self.rst_extension}"
            shutil.copyfile(os.path.join(search_folder, self.rst_file), os.path.join(self.output_bstates_dir, rst_dest_name))
            return rst_dest_name
        elif self.rst_extension:
            possible_hits = [file for file in files if file.endswith(self.rst_extension)]
            if len(possible_hits) > 1:
                possible_hits = sorted(possible_hits, key=lambda file: os.path.getctime(file))
                log.warning(
                    f'Found {possible_hits[-1]} as restart file for iteration {iteration} and walker {walker} based on file creation times. if this is incorrect, provide a file name using flag --rst-file'
                )
            rst_dest_name = f"{iteration:06d}_{walker:06d}.{self.rst_extension}"
            shutil.copyfile(os.path.join(search_folder, possible_hits[-1]), os.path.join(self.output_bstates_dir, rst_dest_name))
            return rst_dest_name
        possible_hits = []
        for test_extension in self.traj_exc_exts:
            possible_hits += [file for file in files if file.endswith(test_extension)]
        if len(possible_hits) < 1:
            for test_extension in self.traj_or_top_exts:
                possible_hits += [file for file in files if file.endswith(test_extension)]
        if len(possible_hits) == 1:
            log.warning(
                f'Found {possible_hits[-1]} as restart file for iteration {iteration} and walker {walker} based on mdtraj extensions. if this is incorrect, provide a file name using flag --rst-file'
            )
        if len(possible_hits) > 1:
            possible_hits = sorted(possible_hits, key=lambda file: os.path.getctime(os.path.join(search_folder, file)))
            log.warning(
                f'Found {possible_hits[-1]} as restart file for iteration {iteration} and walker {walker} based on file creation times. if this is incorrect, provide a file name using flag --rst-file'
            )
        rst_dest_name = f"{iteration:06d}_{walker:06d}.{possible_hits[-1].split('.')[-1]}"
        shutil.copyfile(os.path.join(search_folder, possible_hits[-1]), os.path.join(self.output_bstates_dir, rst_dest_name))
        return rst_dest_name

    def go(self):
        """
        Main public method for running w_reverse. Runs w_succ to find the successful trajectories then iterates over them to copy the trajectories to the output_bstates_dir. Then it iterates over the trajectories again to make the output_bstates_file.
        """
        succ_pairs = self.w_succ()
        # make directory for bstates_reverse if it doesn't already exist
        os.makedirs(self.output_bstates_dir, exist_ok=True)
        # Number of reverse bstates created
        # different totals if the max is less than total succ_pairs to loop
        total_pairs = min(self.max_n_bstates, len(succ_pairs))
        # then for each pair
        rng = np.random.default_rng(self.seed)
        succ_pairs_used = rng.choice(
            succ_pairs, size=total_pairs, p=succ_pairs[:, 2] / np.sum(succ_pairs[:, 2], dtype=float), replace=False
        )
        total_weight = 0
        for idx, succ_pair_used in tqdm(enumerate(succ_pairs_used), total=total_pairs, desc="New bstates"):
            iteration = int(succ_pair_used[0])
            walker = int(succ_pair_used[1])
            weight = float(succ_pair_used[2])
            total_weight += weight
            # check if using HDF5 framework
            if self.h5_framework:
                with tempfile.TemporaryDirectory() as tmpdirname:
                    # Extract the restart data from the .h5 file
                    segment = Segment(n_iter=iteration, seg_id=walker, weight=weight)
                    h5file = WESTIterationFile(self.traj_seg_path.format(n_iter=iteration))
                    h5file.read_restart(segment)
                    restart_writer(tmpdirname, segment)
                    self.copy_traj(tmpdirname, iteration, walker)
            else:
                search_folder = self.traj_seg_path.format(n_iter=iteration, seg_id=walker)
                self.copy_traj(search_folder, iteration, walker)
        # create bstates.txt file
        # fill out the bstates.txt file with name and weight
        # but only use weights if requested, otherwise use equal weights
        # bstates.txt row format: bstate_n | weight | bstate_filename
        with open(os.path.join(self.output_bstates_dir, self.output_bstates_file), "w") as bstates_f:
            for idx, succ_pair_used in enumerate(succ_pairs_used):
                iteration = int(succ_pair_used[0])
                walker = int(succ_pair_used[1])
                weight = float(succ_pair_used[2])
                files = [file for file in os.listdir(self.output_bstates_dir) if file.startswith(f'{iteration:06d}_{walker:06d}')]
                rst_dest_name = ''
                if len(files) > 1:
                    if f'{iteration:06d}_{walker:06d}.{self.rst_extension}' in files:
                        rst_dest_name = f'{iteration:06d}_{walker:06d}.{self.rst_extension}'
                    else:
                        files = sorted(files, key=lambda file: os.path.getctime(file))
                        rst_dest_name = files[-1]
                        log.warning(
                            f'Muliple files starting with {iteration:06d}_{walker:06d} are in the output directory. Using {rst_dest_name}. if this is incorrect, provide a file name using flag --rst-file so that the correct extension can be used'
                        )
                elif len(files) == 1:
                    rst_dest_name = files[0]
                else:
                    log.error(
                        f'A restart file starting with {iteration:06d}_{walker:06d} should be present in the output directory but is not!!!'
                    )
                if self.use_weights:
                    bstates_f.write(f"{idx} {weight / total_weight:.3e} {rst_dest_name}\n")
                else:
                    bstates_f.write(f"{idx} {1 / total_pairs:.3e)} {rst_dest_name}\n")


def entry_point():
    W_Reverse().main()


if __name__ == "__main__":
    entry_point()
