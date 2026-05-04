import logging
import h5py
from tqdm.auto import tqdm
import os
import shutil
import argparse
from io import BytesIO
import tarfile
from westpa.core.h5io import WESTIterationFile
from westpa.core.h5io import safe_extract
from westpa.core._rc import WESTRC
from westpa.tools import (
    WESTTool,
    WESTDataReader,
    IterRangeSelection,
)
log = logging.getLogger('w_reverse')


class W_Reverse:
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

    def __init__(
        self,
        h5="west.h5",
        first_iter=1,
        last_iter=None,
        config_file="west.cfg",
        max_n_bstates=10000,
        rst_file='seg.ncrst',
        output_bstates_dir="bstates_reverse",
        output_bstates_file="bstates.txt",
        use_weights=True,
        temp_dir="temp_dir",
    ):
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
        Name of the temporary directory that will be created
        """
        # Parse config file
        westrc = WESTRC()
        westrc.read_config(config_file)
        self.config_file = config_file
        config = westrc.config
        # Read the west.h5 file
        self.h5 = h5py.File(h5, mode="r")
        self.first_iter = int(first_iter)
        # default to last
        if last_iter is not None:
            self.last_iter = int(last_iter)
        elif last_iter is None:
            self.last_iter = self.h5.attrs["west_current_iteration"] - 1
        # Look at the data_refs from the config file
        data_refs_dic = config['west']['data']['data_refs']
        # Default to not using HDF5 framework
        self.h5_framework = False
        if 'iteration' in data_refs_dic.keys():
            traj_seg_path_list = data_refs_dic['iteration'].split('/')[1:-1]
            self.h5_framework = True
        else:
            traj_seg_path_list = data_refs_dic['segment'].split('/')[1:]
        self.traj_segs_path = '/'.join(traj_seg_path_list)
        self.max_n_bstates = int(max_n_bstates)
        self.rst_file = str(rst_file)
        # Get the restart file extention being used
        self.rst_extension = self.rst_file.split('.')[-1]
        self.output_bstates_dir = str(output_bstates_dir)
        self.temp_dir = str(temp_dir)
        self.output_bstates_file = str(output_bstates_file)
        self.use_weights = use_weights

    def w_succ(self):
        """
        Find and return all successfully recycled (iter, seg) pairs.

        Returns
        -------
        succ : list of tuples (iter,wlk,weight)
        """
        succ = []

        for iter in tqdm(range(self.first_iter, self.last_iter + 1), desc="Running w_succ"):
            # if the new_weights group exists in the h5 file
            if f"iterations/iter_{iter:08d}/new_weights" in self.h5:
                prev_segs = self.h5[f"iterations/iter_{iter:08d}/new_weights/index"]["prev_seg_id"]
                recycled_weights = self.h5[f"iterations/iter_{iter:08d}/new_weights/index"]["weight"]
                # append the previous iter and previous seg id recycled and the weight
                for i in range(len(prev_segs)):
                    succ.append((iter - 1, prev_segs[i], recycled_weights[i]))
        return succ

    @staticmethod
    def create_dir(directory):
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Directory '{directory}' created.")
        else:
            print(f"Directory '{directory}' already exists.")

    def w_reverse(self):
        """
        Main public method for running w_reverse.
        """
        # Start a new instance of WESTRC to parse the config file
        westrc = WESTRC()
        westrc.read_config(self.config_file)
        config = westrc.config
        # Get the data from the data_refs
        data_refs_dic = config['west']['data']['data_refs']
        # first generate list of successful events / iter,seg pairs and weights
        succ_pairs = self.w_succ()
        # Data I was used for testing
        # succ_pairs = [(73, 130, 5.991585103556223e-13), (74, 132, 7.489481379445279e-14), (74, 150, 7.489481379445279e-14)]
        # succ_pairs = [(73, 130, 5.991585103556223e-13)]
        # make directory for bstates_reverse if it doesn't already exist
        self.create_dir(self.output_bstates_dir)
        # make directory for tmporary file storage
        self.create_dir(self.temp_dir)
        # create bstates.txt file
        with open(f"{self.output_bstates_dir}/{self.output_bstates_file}", "w") as bstates_f:

            # Number of reverse bstates created
            n_bstates = 0
            # different totals if the max is less than total succ_pairs to loop
            if self.max_n_bstates < len(succ_pairs):
                total_pairs = self.max_n_bstates
            else:
                total_pairs = len(succ_pairs)

            # then for each pair
            for idx, (it, wlk, weight) in tqdm(enumerate(succ_pairs), total=total_pairs, desc="New bstates"):
                # Make sure you are not over the maximum bstates
                if n_bstates < self.max_n_bstates:
                    # Assign new bstate restart file name
                    rst_dest_name = f"{it:06d}_{wlk:06d}.{self.rst_extension}"
                    # check if using HDF5 framework
                    if self.h5_framework:
                        # Find how the .h5 files for each iteration are named
                        traj_seg_file_name = data_refs_dic["iteration"].split('/')[-1]
                        # Make a path to this iterations .h5 file
                        traj_seg = f'{self.traj_segs_path}/{traj_seg_file_name}'
                        # Extracct the restart data from the .h5 file
                        h5file = WESTIterationFile(traj_seg.format(n_iter=it))
                        restart_data = h5file.read_data('/restart/%d_%d' % (it, wlk), 'data')
                        try:
                            if restart_data is None:
                                raise ValueError('restart data is not present')
                            # Extract all restart files into a temporary directory
                            with BytesIO(restart_data[:-1]) as d:
                                with tarfile.open(fileobj=d, mode='r:gz') as t:
                                    safe_extract(t, path=self.temp_dir)
                        except ValueError as e:
                            log.warning(f'could not write HDF5 Framework restart data for iteration {it} walker {wlk}: {e}')
                            if it == 1:
                                log.warning(
                                    f'In iteration 1. Assuming this is a start state and proceeding to skip reading restart from per-iteration HDF5 file for iteration {it} walker {wlk}'
                                )
                        except Exception as e:
                            log.warning(f'could not write HDF5 Framework restart data for iteration {it} walker {wlk}: {e}')
                        # Look at all files in the temp directory
                        temp_dir_contents = os.listdir(self.temp_dir)
                        extention_not_found = True
                        for temp_file in temp_dir_contents:
                            # Move and rename only the restart file with the user specified extension to the bstate directory
                            if temp_file.split('.')[-1] == self.rst_extension:
                                extention_not_found = False
                                shutil.move(
                                    f"{self.temp_dir}/{temp_file}",
                                    f"{self.output_bstates_dir}/{it:06d}_{wlk:06d}.{self.rst_extension}",
                                )
                                break
                        if extention_not_found:
                            log.warning(
                                f'File with extension {self.rst_extension} is not present in the restart data of {traj_seg.format(n_iter=it)}'
                            )
                    else:
                        # find the corresponding restart file
                        seg_path = (
                            self.traj_segs_path.replace('segment.n_iter', 'n_iter')
                            .replace('segment.seg_id', 'seg_id')
                            .format(n_iter=it, seg_id=wlk)
                        )
                        os.listdir(seg_path)
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
                    n_bstates += 1
                else:
                    break
            # Remove temporary directory
            shutil.rmtree(self.temp_dir)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="w_reverse: a tool for generating bstates for a subsequent "
        + "steady-state WE simulation in the opposite direction."
    )
    parser.add_argument(
        "-W", "-w", "--west", "--west-data", "-h5", "--h5file", dest="h5", type=str, default="west.h5", help="Path to west.h5 file"
    )
    parser.add_argument(
        "--first-iter", "-fi", dest="first_iter", type=int, default=1, help="First iteration to consider (default: 1)"
    )
    parser.add_argument(
        "--last-iter",
        "-li",
        dest="last_iter",
        type=int,
        default=None,
        help="Last iteration to consider (default: last recorded iteration in west.h5)",
    )
    parser.add_argument("--config-file", dest="config_file", type=str, default="west.cfg", help="Path to the config file")
    parser.add_argument(
        "--max-n-bstates",
        dest="max_n_bstates",
        type=int,
        default=10000,
        help="Max number of bstates to copy over. Adjust this if you prefer "
        + "a subset of the first bstates found. Default max of 10000.",
    )
    parser.add_argument("--rst-file", dest="rst_file", type=str, default="seg.ncrst", help="Path to the Restart File")
    parser.add_argument(
        "--output-bstates-dir",
        "-obd",
        dest="output_bstates_dir",
        type=str,
        default="bstates_reverse",
        help="Output directory for the bstates and output_bstates_file",
    )
    parser.add_argument(
        "--output-bstates-file",
        "-obf",
        dest="output_bstates_file",
        type=str,
        default="bstates.txt",
        help="Name of the output bstates file",
    )
    # TODO: may need to be adjusted to store False when included
    parser.add_argument(
        "--no-weights",
        "-nw",
        dest="use_weights",
        action="store_false",
        help="Don't include the recycled event weight when making the bstates.txt file",
    )
    parser.add_argument(
        "--temp-dir",
        "-td",
        dest="temp_dir",
        type=str,
        default="temp_dir",
        help="Directory that files will temporarily be stored in while w_reverse is running",
    )
    return parser.parse_args()


def entry_point():
    args = parse_arguments()
    reverse = W_Reverse(**vars(args))
    reverse.w_reverse()


if __name__ == "__main__":
    entry_point()
