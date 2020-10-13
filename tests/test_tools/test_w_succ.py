import os
import filecmp
import argparse

from westpa.cli.core.w_succ import entry_point
from unittest import mock


class Test_W_Succ:
    '''
    Class to test that w_succ is working successfully to find and output recycled trajectories.
    '''

    def test_run_w_succ(self, ref_initialized):

        w_succ_ref = os.path.join(self.REFERENCE_PATH, 'w_succ_ref')
        temp_w_succ_output_file = open('w_succ_temp', 'w')

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity=0,
                rcfile=self.cfg_filepath,
                output_file=temp_w_succ_output_file,
                west_h5name=self.h5_filepath,
                anal_h5name=self.h5_filepath,
            ),
        ):
            entry_point()

        temp_w_succ_output_file.close()

        assert os.path.isfile('w_succ_temp'), 'The output file was not generated.'
        assert filecmp.cmp(w_succ_ref, 'w_succ_temp'), 'The output file was not generated correctly'

        os.remove('w_succ_temp')
