
import nose
import nose.tools
from nose.tools import with_setup
from nose.plugins.skip import SkipTest

from .common import *
from w_trace import WTraceTool
import h5py


class Test_W_Trace_Args(CommonToolTest):
    '''Class to test w_trace works with different positional arguments (of the form n_iter:seg_id [n_iter:seg_id])
    This class tests that a) w_trace works with different numbers of positional arguments, and b) that the output file is set up
    correctly AND the data from subsequent traces is appended (in other words, if there is an existing output file, it
    is not overwritten each time trace is called)'''

    arg_combos = [['20:0'], ['20:1', '20:2']]

    test_name = 'W_TRACE'

    def test_args(self):
        '''Testing arg combos: w_trace runs as expected'''

        self.args_so_far = []
        self.endpoints = []
        self.outfile = self.mktemp(prefix='trace')
        outarg = '-o={}'.format(self.outfile)

        for args in self.arg_combos:
            self.w = WTraceTool()
            
            self.args_so_far.extend(args)
            args.append(outarg)
            self.w.make_parser_and_process(args=args)
            self.endpoints.extend(self.w.endpoints)
            args.remove(outarg)

            test_outcome, err = self.check_runs_with_args()

            yield self.check_args, test_outcome, err, args

            del self.w

        yield self.check_args, test_outcome, err, args

    def check_output(self):

        with h5py.File(self.outfile) as f:

            assert 'trajectories' in list(f.keys()), "'trajectories' group not in output file"

            traj_group = f['trajectories']

            #Expected trace groups - successive runs of w_trace with different traces should add
            # newly traced groups to output hdf5 file
            expected_groups = sorted([self.w.output_pattern % (n_iter,seg_id) for n_iter, seg_id in self.endpoints])

            assert list(traj_group.keys()) == expected_groups, "H5 groups ({}) are not as expected ({})".format(list(traj_group.keys()), expected_groups)