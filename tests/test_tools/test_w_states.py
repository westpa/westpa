import argparse

from westpa.cli.core.w_states import entry_point
from unittest import mock
from pytest import fixture


class Test_W_States:
    def test_run_w_states(self, ref_initialized):
        '''Tests running w_states on a sample h5 file'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(mode='show', verbosity=0, rcfile=self.cfg_filepath, bstate_file=None),
        ):
            entry_point()

        out, err = self.capfd.readouterr()

        state_string = out.split('\n')[-2]
        known_state_string = 'initial                1    None        # state_id=0    pcoord=[8.0]'

        assert state_string == known_state_string

    @fixture(autouse=True)
    def capfd(self, capfd):
        self.capfd = capfd
