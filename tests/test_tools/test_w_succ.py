import os
import filecmp
from westpa.cli.core.w_succ import WSucc
import westpa
from common import MockArgs


class Test_W_Succ:
    '''
    Class to test that w_succ is working successfully to find and output recycled trajectories.
    '''

    def test_run_w_succ(self, ref_50iter):
        temp_w_succ_output_file = open('w_succ_temp.txt', 'w')

        args = MockArgs(
            rcfile=self.cfg_filepath,
            verbosity='debug',
            output_file=temp_w_succ_output_file,
            west_h5name=self.h5_filepath,
            anal_h5name='analysis.h5',
        )

        westpa.rc.process_args(args, config_required=False)

        w_succ = WSucc()
        w_succ.process_args(args)
        w_succ.output_file = temp_w_succ_output_file

        w_succ.find_successful_trajs()

        temp_w_succ_output_file.close()
        assert os.path.isfile('./w_succ_temp.txt'), 'The output file was not generated.'
        assert filecmp.cmp('w_succ_ref.txt', 'w_succ_temp.txt'), 'The output file was not generated correctly'

        os.remove('w_succ_temp.txt')
