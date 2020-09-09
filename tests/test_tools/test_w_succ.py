import os
import filecmp
import unittest
import argparse

from westpa.cli.core.w_succ import entry_point
from unittest import mock


class Test_W_Succ(unittest.TestCase):
    '''
    Class to test that w_succ is working successfully to find and output recycled trajectories.
    '''

    test_name = 'W_SUCC'

    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_succ')

    def setUp(self):

        self.starting_path = os.getcwd()
        self.odld_path = os.path.dirname(__file__) + '/refs'
        os.chdir(self.odld_path)

    def test_run_w_succ(self):

        temp_w_succ_output_file = open('w_succ_temp', 'w+')

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile='west.cfg',
                output_file=temp_w_succ_output_file,
                west_h5name='west.h5',
                anal_h5name='analysis.h5',
            ),
        ):
            entry_point()

        temp_w_succ_output_file.close()

        assert os.path.isfile('w_succ_temp'), 'The output file was not generated.'
        self.assertTrue(filecmp.cmp('w_succ_ref', 'w_succ_temp'), 'The output file was not generated correctly')

    def tearDown(self):

        os.remove(self.odld_path + '/w_succ_temp')
        os.chdir(self.starting_path)
