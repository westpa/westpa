import argparse

from westpa import rc

from westpa.cli.core.w_truncate import entry_point
from unittest import mock
import h5py


class Test_W_Truncate:
    test_name = 'W_TRUNCATE'

    def test_run_w_truncate(self, w_truncate_fixture):
        '''Tests running w_truncate on a 3 iteration h5 file'''

        _hfile = h5py.File(self.h5_filepath, mode='r')
        assert 'iter_00000003' in _hfile['/iterations'].keys()
        _hfile.close()

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(verbosity='debug', rcfile=self.cfg_filepath, n_iter=2),
        ):
            entry_point()

        _hfile = h5py.File(self.h5_filepath, mode='r')
        assert 'iter_00000003' not in _hfile['/iterations'].keys()
        assert 'iter_00000002' not in _hfile['/iterations'].keys()
        _hfile.close()

    def tearDown(self):

        rc._sim_manager = None
        rc._system = None
        rc._data_manager = None
        rc._we_driver = None
        rc._propagator = None
