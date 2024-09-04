import h5py
import argparse
import pytest
from unittest import mock

from westpa.cli.core.w_truncate import entry_point


class Test_W_Truncate:
    def test_run_w_truncate(self, ref_3iter):
        '''Tests running w_truncate on a 3 iteration h5 file'''

        with h5py.File(self.h5_filepath, mode='r') as _hfile:
            assert 'iter_00000003' in list(_hfile['/iterations'].keys())

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(verbosity='debug', rcfile=self.cfg_filepath, n_iter=2, we_h5filename=None),
        ):
            entry_point()

        with h5py.File(self.h5_filepath, mode='r') as _hfile:
            assert 'iter_00000003' not in list(_hfile['/iterations'].keys())
            assert 'iter_00000002' not in list(_hfile['/iterations'].keys())
            assert _hfile.attrs['west_current_iteration'] == 1
        _hfile.close()

    def test_run_w_truncate_default(self, ref_3iter):
        '''Tests running w_truncate on a 3 iteration h5 file, but running with n_iter=0'''

        with h5py.File(self.h5_filepath, mode='r') as _hfile:
            assert 'iter_00000003' in list(_hfile['/iterations'].keys())

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(verbosity='debug', rcfile=self.cfg_filepath, n_iter=0, we_h5filename=None),
        ):
            entry_point()

        with h5py.File(self.h5_filepath, mode='r') as _hfile:
            assert 'iter_00000003' not in list(_hfile['/iterations'].keys())
            assert 'iter_00000002' in list(_hfile['/iterations'].keys())
            assert _hfile.attrs['west_current_iteration'] == 2

    @pytest.mark.xfail(raises=SystemExit)
    def test_run_w_truncate_fail(self, ref_3iter):
        '''Tests running w_truncate on a 3 iteration h5 file, but running with n_iter=0'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(verbosity='debug', rcfile=self.cfg_filepath, n_iter=5, we_h5filename=None),
        ):
            entry_point()
