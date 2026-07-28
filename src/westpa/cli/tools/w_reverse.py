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
from numpy import flatnonzero, random, array, sum
from westpa.core.trajectory import find_top_traj_file, mdtraj_supported_extensions

log = logging.getLogger('w_reverse')


class W_Reverse(WESTTool):
    """
    w_reverse: a tool for taking a WE simulation facilitated
    through WESTPA with successful recycling events and generating
    a new directory with the recycled restart files, which can then
    serve as the bstates for a subsequent WE simulation in the opposite
    direction, i.e. starting from the successfully recycled events and
    then going back to the original starting point; overall reversed.

    TODO:
        * option to use w_assign/assign.h5 output for successful bstate selection
        * adapt for WE simulation traj_segs that used the hdf5 framework
        * integrate into WESTPA, pull instance attributes and args from west.cfg
        * change printing to west logging
    """

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
        # TODO: may need to be adjusted to store False when included
        rgroup.add_argument(
            "--no-weights",
            "-nw",
            dest="use_weights",
            action="store_false",
            help="Don't include the recycled event weight when making the bstates.txt file",
        )

    def process_args(self, args):
        """
        Parameters
        ----------
        h5 : str
            Path to west.h5 file
        first_iter : int
            By default start at iteration 1.
        last_iter : int
            Last iteration data to include, default is the last recorded iteration in the west.h5 file.
        config_file : str
            Name of the configuration file
            max_n_bstates : int
            Max number of bstates to copy over. Adjust this if you prefer only the first
            n bstates found, default 10,000.
        rst_file : str
            Name of the restart file within each traj_segs/ subdirectory.
        output_bstates_dir : str
            Output directory for the bstates and output_bstates_file.
            Default './bstates_reverse'.'
        output_bstates_file : str
            Name of the output bstates file, default 'bstates.txt'.
        use_weights : bool
            By default, include the recycled event weight when making the bstates.txt file.
            temp_dir : str
        """
        self.data_reader.process_args(args)
        self.config_required = True
        self.config_file = args.config_file
        self.westrc.read_config(self.config_file)
        self.config = self.westrc.config
        # Read the west.h5 file
        self.h5 = WESTPAH5File(args.we_h5filename, 'r')
        self.first_iter = int(args.first_iter)
        # default to last
        if args.last_iter is not None:
            self.last_iter = int(args.last_iter)
        elif args.last_iter is None:
            self.last_iter = self.h5.attrs['west_current_iteration'] - 1
        # Look at the data_refs from the config file
        self.data_refs_dic = self.config['west']['data']['data_refs']
        # Default to not using HDF5 framework
        self.h5_framework = False
        starts_with_slash = False
        if 'iteration' in self.data_refs_dic.keys():
            traj_seg_file_name = self.data_refs_dic['iteration'].split('/')[-1]
            if self.data_refs_dic['iteration'][0] == '$':
                traj_seg_path_list = self.data_refs_dic['iteration'].split('/')[1:-1]
            else:
                traj_seg_path_list = self.data_refs_dic['iteration'].split('/')[:-1]
            if self.data_refs_dic['iteration'][0] == '/':
                starts_with_slash = True
            self.h5_framework = True
        else:
            if self.data_refs_dic['segment'][0] == '$':
                traj_seg_path_list = self.data_refs_dic['segment'].split('/')[1:]
            else:
                traj_seg_path_list = self.data_refs_dic['segment'].split('/')
            if self.data_refs_dic['segment'][0] == '/':
                starts_with_slash = True
        self.traj_segs_path = '/'.join(traj_seg_path_list)
        if starts_with_slash:
            self.traj_segs_path = f'/{self.traj_segs_path}'
        if self.h5_framework:
            self.traj_seg = f'{self.traj_segs_path}/{traj_seg_file_name}'
        self.max_n_bstates = int(args.max_n_bstates)
        if args.rst_file:
            self.rst_file = str(args.rst_file).lower()
            # Get the restart file extension being used
            self.rst_extension = self.rst_file.split('.')[-1]
        else:
            self.rst_file = None
            self.rst_extension = None
        self.traj_exc_exts = []
        self.traj_or_top_exts = []
        for i in self.traj_exts:
            if i in self.top_exts:
                self.traj_or_top_exts.append(i)
            else:
                self.traj_exc_exts.append(i)
        self.output_bstates_dir = str(args.output_bstates_dir)
        if os.path.isdir(self.output_bstates_dir):
            shutil.rmtree(self.output_bstates_dir)
        self.output_bstates_file = str(args.output_bstates_file)
        self.use_weights = args.use_weights

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
            endpoint_type = self.h5['iterations'][iteration]['seg_index']['endpoint_type']
            indices = flatnonzero(endpoint_type == Segment.SEG_ENDPOINT_RECYCLED)
            temp_array = [
                [iteration_index + 1, index, self.h5['iterations'][iteration]['seg_index']['weight'][index]] for index in indices
            ]
            succ += temp_array
        return array(succ)

    def go(self):
        """
        Main public method for running w_reverse.
        """
        succ_pairs = self.w_succ()
        # Data I was used for testing
        # succ_pairs = [(73, 130, 5.991585103556223e-13), (74, 132, 7.489481379445279e-14), (74, 150, 7.489481379445279e-14)]
        # succ_pairs = [(73, 130, 5.991585103556223e-13)]
        # make directory for bstates_reverse if it doesn't already exist
        os.makedirs(self.output_bstates_dir)
        # create bstates.txt file
        with open(f"{self.output_bstates_dir}/{self.output_bstates_file}", "w") as bstates_f:

            # Number of reverse bstates created
            # different totals if the max is less than total succ_pairs to loop
            if self.max_n_bstates < len(succ_pairs):
                total_pairs = self.max_n_bstates
            else:
                total_pairs = len(succ_pairs)
            # then for each pair
            rng = random.default_rng(12345)
            indices = rng.choice(
                len(succ_pairs), size=total_pairs, p=succ_pairs[:, 2] / sum(succ_pairs[:, 2], dtype=float), replace=False
            )
            for idx, index in tqdm(enumerate(indices), total=total_pairs, desc="New bstates"):
                iteration = int(succ_pairs[index][0])
                walker = int(succ_pairs[index][1])
                weight = float(succ_pairs[index][2])
                # check if using HDF5 framework
                if self.h5_framework:
                    with tempfile.TemporaryDirectory() as tmpdirname:
                        # Extracct the restart data from the .h5 file
                        h5file = WESTIterationFile(self.traj_seg.format(n_iter=iteration))
                        restart_data = h5file.read_data('/restart/%d_%d' % (iteration, walker), 'data')
                        segment = Segment(n_iter=iteration, seg_id=walker, weight=weight)
                        segment.data['iterh5/restart'] = restart_data
                        restart_writer(tmpdirname, segment)
                        # Look at all files in the temp directory
                        if self.rst_file:
                            rst_dest_name = f"{iteration:06d}_{walker:06d}.{self.rst_extension}"
                            temp_dir_contents = os.listdir(tmpdirname)
                            file_not_found = True
                            for temp_file in temp_dir_contents:
                                # Move and rename only the restart file with the user specified extension to the bstate directory
                                if temp_file.lower() == self.rst_file:
                                    file_not_found = False
                                    shutil.move(
                                        f"{tmpdirname}/{temp_file}",
                                        f"{self.output_bstates_dir}/{iteration:06d}_{walker:06d}.{self.rst_extension}",
                                    )
                                    break
                            if file_not_found:
                                log.warning(
                                    f"File {self.rst_file} is not present in the restart data of {self.traj_seg.format(n_iter=iteration)}"
                                )
                        else:
                            _, traj_file = find_top_traj_file(tmpdirname, [], self.traj_exc_exts)
                            if not traj_file:
                                _, traj_file = find_top_traj_file(tmpdirname, [], self.traj_or_top_exts)
                            extension = traj_file.split('/')[-1].split('.')[-1].lower()
                            rst_dest_name = f"{iteration:06d}_{walker:06d}.{extension}"
                            shutil.move(
                                traj_file,
                                f"{self.output_bstates_dir}/{rst_dest_name}",
                            )
                else:
                    if not self.rst_file:
                        log.warning(
                            "The flag --rst-file that defines the restart file name must be used if the HDF5 frame work is not used!!"
                        )
                    rst_dest_name = f"{iteration:06d}_{walker:06d}.{self.rst_extension}"
                    # find the corresponding restart file
                    seg_path = (
                        self.traj_segs_path.replace('segment.n_iter', 'n_iter')
                        .replace('segment.seg_id', 'seg_id')
                        .format(n_iter=iteration, seg_id=walker)
                    )
                    # os.listdir(seg_path)
                    rst_file_path = f'{seg_path}/{self.rst_file}'
                    # if bstate file exists, skip
                    if os.path.exists(f"{self.output_bstates_dir}/{rst_dest_name}"):
                        continue
                    shutil.copyfile(rst_file_path, f"{self.output_bstates_dir}/{rst_dest_name}")
                # fill out the bstates.txt file with name and weight
                # but only use weights if requested, otherwise use equal weights
                # bstates.txt row format: bstate_n | weight | bstate_filename
                if self.use_weights:
                    bstates_f.write(f"{idx} {weight:.3e} {rst_dest_name}\n")
                else:
                    bstates_f.write(f"{idx} 1 {rst_dest_name}\n")


def entry_point():
    W_Reverse().main()


if __name__ == "__main__":
    entry_point()
