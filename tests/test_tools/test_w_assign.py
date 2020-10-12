import unittest
import argparse
import os
import shutil

from h5diff import H5Diff

from westpa.cli.tools.w_assign import entry_point
from unittest import mock


class Test_W_Assign(unittest.TestCase):

    test_name = 'W_ASSIGN'

    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_assign')

    def setUp(self):

        self.starting_path = os.getcwd()

        self.ref_path = os.path.dirname(__file__) + '/ref'
        os.environ['WEST_SIM_ROOT'] = self.ref_path

        os.chdir(self.ref_path)

    def test_run_w_assign(self):

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile='west.cfg',
                max_queue_length=None,
                we_h5filename='west.h5',
                construct_dataset=None,
                dsspecs=None,
                output='assign.h5',
                subsample=None,
                config_from_file=True,
                scheme='TEST',
            ),
        ):

            entry_point()

        diff = H5Diff(self.ref_path + '/assign_ref.h5', self.ref_path + '/ANALYSIS/TEST/assign.h5')
        diff.check()

    def tearDown(self):

        os.environ['WEST_SIM_ROOT'] = ''
        shutil.rmtree('ANALYSIS')
        os.chdir(self.starting_path)
