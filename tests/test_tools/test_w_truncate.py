import argparse

from westpa.cli.core.w_truncate import entry_point
from unittest import mock
import h5py

import os


class Test_W_Truncate:
    def test_run_w_truncate(self, ref_3iter):
        '''Tests running w_truncate on a 3 iteration h5 file'''

        _hfile = h5py.File(self.h5_filepath, mode='r')
        assert 'iter_00000003' in list(_hfile['/iterations'].keys())
        _hfile.close()
        del _hfile

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(verbosity='debug', rcfile=self.cfg_filepath, n_iter=2),
        ):
            entry_point()

        print(os.getcwd())
        _hfile = h5py.File(self.h5_filepath, mode='r')
        print(_hfile['/iterations'].keys())
        assert 'iter_00000003' not in list(_hfile['/iterations'].keys())
        assert 'iter_00000002' not in list(_hfile['/iterations'].keys())
        _hfile.close()
