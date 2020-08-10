import unittest
import argparse
import os

import westpa

from westpa.cli.core.w_truncate import entry_point
from unittest import mock
import h5py

from shutil import copyfile


class Test_W_Truncate(unittest.TestCase):
    test_name = 'W_TRUNCATE'

    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_truncate')

    def setUp(self):

        self.starting_path = os.getcwd()

        self.odld_path = os.path.dirname(__file__) + '/ref'

    def test_run_w_truncate(self):
        '''Tests running w_truncate on a 3 iteration h5 file'''

        self.starting_path = os.getcwd()

        odld_path = os.path.dirname(__file__) + '/ref'

        os.chdir(odld_path)
        copyfile('west_3iter.h5', 'west.h5')

        _hfile = h5py.File(self.odld_path + '/west.h5', mode='r')
        assert 'iter_00000003' in _hfile['/iterations'].keys()
        _hfile.close()

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(verbosity='debug', rcfile='west.cfg', n_iter=2),
        ):
            entry_point()

        _hfile = h5py.File(self.odld_path + '/west.h5', mode='r')
        assert 'iter_00000003' not in _hfile['/iterations'].keys()
        assert 'iter_00000002' not in _hfile['/iterations'].keys()
        _hfile.close()

    def tearDown(self):

        westpa.rc._sim_manager = None
        westpa.rc._system = None
        westpa.rc._data_manager = None
        westpa.rc._we_driver = None
        westpa.rc._propagator = None
        os.environ['WEST_SIM_ROOT'] = ''

        os.remove(self.odld_path + '/west.h5')
        os.chdir(self.starting_path)
