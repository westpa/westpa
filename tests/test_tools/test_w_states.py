import unittest
import argparse
import os

import westpa

from westpa.cli.core.w_states import entry_point
from unittest import mock
from pytest import fixture

from shutil import copyfile


class Test_W_States(unittest.TestCase):
    test_name = 'W_STATES'

    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_states')

    def setUp(self):

        self.starting_path = os.getcwd()

        self.odld_path = os.path.dirname(__file__) + '/ref'

    def test_run_w_states(self):
        '''Tests running w_states on a sample h5 file'''

        os.chdir(self.odld_path)
        copyfile('west_ref.h5', 'west.h5')

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(mode='show', verbosity=0, rcfile='west.cfg', bstate_file=None),
        ):
            entry_point()

        out, err = self.capfd.readouterr()

        state_string = out.split('\n')[-2]
        known_state_string = 'initial                1    None        # state_id=0    pcoord=[8.0]'

        assert state_string == known_state_string

    @fixture(autouse=True)
    def capfd(self, capfd):
        self.capfd = capfd

    def tearDown(self):

        westpa.rc._sim_manager = None
        westpa.rc._system = None
        westpa.rc._data_manager = None
        westpa.rc._we_driver = None
        westpa.rc._propagator = None
        os.environ['WEST_SIM_ROOT'] = ''

        os.remove(self.odld_path + '/west.h5')
        os.chdir(self.starting_path)
