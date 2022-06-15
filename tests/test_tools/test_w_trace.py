import unittest

import h5py
import pytest
import tempfile
import argparse
import os

from common import CommonToolTest
from westpa.cli.tools.w_trace import WTraceTool, entry_point
from h5diff import H5Diff
from filecmp import cmp
from unittest import mock


class Test_NEW_W_Trace:
    '''Class to test w_trace by running through the entry_point and comparing the output text files with reference files.'''

    def test_trace(self, ref_50iter):
        arg_combos = [['20:0'], ['20:1'], ['20:2'], ['20:1', "20:2"]]
        output_combos = ['traj_20_0_trace.txt', 'traj_20_1_trace.txt', 'traj_20_2_trace.txt', 'traj_20_2_trace.txt']
        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')

        for arg, output_txt in zip(arg_combos, output_combos):
            self.outfile = tempfile.TemporaryFile(suffix='.h5', prefix='trace')
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

            self.outfile.close()

            # Compare text file output
            assert cmp(
                os.path.join(test_dir, output_txt), os.path.join(ref_dir, output_txt), shallow=False
            ), f'Output file {output_txt} is not the same as reference file.'


@pytest.mark.skip(reason='doesn\'t actually work')
class Test_W_Trace_Args(unittest.TestCase, CommonToolTest):
    '''Class to test w_trace works with different positional arguments (of the form n_iter:seg_id [n_iter:seg_id])
    This class tests that a) w_trace works with different numbers of positional arguments, b) that the output file is set up
    correctly AND the data from subsequent traces is appended (in other words, if there is an existing output file, it
    is not overwritten each time trace is called), and c) that the content of the h5 files generated match the content of previously
    generated reference files.'''

    arg_combos = [['20:0'], ['20:1', '20:2']]
    ref_files = ['refs/trajs_ref_20_0.h5', 'refs/trajs_ref_full.h5']

    test_name = 'W_TRACE'

    def test_args(self):
        '''Testing arg combos: w_trace runs as expected'''
        self.args_so_far = []
        self.endpoints = []
        self.outfile = self.mktemp(prefix='trace')
        outarg = '-o={}'.format(self.outfile)

        for idx, args in enumerate(self.arg_combos):
            self.w = WTraceTool()

            self.args_so_far.extend(args)
            args.append(outarg)
            self.w.make_parser_and_process(args=args)
            self.endpoints.extend(self.w.endpoints)
            args.remove(outarg)

            test_outcome, err = self.check_runs_with_args(idx)

            return self.check_args, test_outcome, err, args

            del self.w

        return self.check_args, test_outcome, err, args

    def check_output(self, arg_idx):

        with h5py.File(self.outfile) as f:

            assert 'trajectories' in list(f.keys()), "'trajectories' group not in output file"

            traj_group = f['trajectories']

            # Expected trace groups - successive runs of w_trace with different traces should add
            # newly traced groups to output hdf5 file
            expected_groups = sorted([self.w.output_pattern % (n_iter, seg_id) for n_iter, seg_id in self.endpoints])

            assert list(traj_group.keys()) == expected_groups, "H5 groups ({}) are not as expected ({})".format(
                list(traj_group.keys()), expected_groups
            )

        diff = H5Diff(self.ref_files[arg_idx], self.outfile)
        diff.check()
