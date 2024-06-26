import h5py
import pytest
import tempfile
import argparse
import os

from westpa.cli.tools.w_trace import entry_point
from h5diff import H5Diff
from filecmp import cmp
from unittest import mock


class Test_NEW_W_Trace:
    '''Class to test w_trace works with different positional arguments (of the form n_iter:seg_id [n_iter:seg_id])
    This class tests that a) w_trace works with different numbers of positional arguments, b) that the output file is set up
    correctly AND the data from subsequent traces is appended (in other words, if there is an existing output file, it
    is not overwritten each time trace is called), and c) that the content of the h5 files and text files generated match the
    content of previously generated reference files.'''

    @pytest.mark.parametrize(
        ['arg', 'output_txt', 'ref_file'],
        [
            [['20:0'], 'traj_20_0_trace.txt', 'trajs_ref_20_0.h5'],
            [['20:1'], 'traj_20_1_trace.txt', None],
            [['20:2'], 'traj_20_2_trace.txt', None],
            [['20:1', "20:2"], 'traj_20_2_trace.txt', 'trajs_ref_full.h5'],
        ],
    )
    def test_trace(self, ref_50iter, arg, output_txt, ref_file):
        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')

        self.outfile = tempfile.TemporaryFile(suffix='.h5', prefix='trace')
        self.endpoints = [[int(z.split(':')[0]), int(z.split(':')[1])] for z in arg]
        self.ref_file = ref_file
        self.output_pattern = 'traj_%d_%d'

        test_dir = os.getcwd()
        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                rcfile=self.cfg_filepath,
                we_h5filename=self.h5_filepath,
                endpoints=arg,
                output=self.outfile,
                verbosity='debug',
                output_pattern='traj_%d_%d',
                datasets=None,
            ),
        ):
            entry_point()

        if ref_file:
            self.check_output()

        self.outfile.close()

        # Compare text file output
        assert cmp(
            os.path.join(test_dir, output_txt), os.path.join(ref_dir, output_txt), shallow=False
        ), f'Output file {output_txt} is not the same as reference file.'

    def check_output(self):
        with h5py.File(self.outfile) as f:
            assert 'trajectories' in list(f.keys()), "'trajectories' group not in output file"

            traj_group = f['trajectories']

            # Expected trace groups - successive runs of w_trace with different traces should add
            # newly traced groups to output hdf5 file
            expected_groups = sorted([self.output_pattern % (n_iter, seg_id) for n_iter, seg_id in self.endpoints])

            assert list(traj_group.keys()) == expected_groups, "H5 groups ({}) are not as expected ({})".format(
                list(traj_group.keys()), expected_groups
            )

        diff = H5Diff(self.ref_file, self.outfile)
        diff.check()
