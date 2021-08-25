from pytest import fixture
import westpa.cli.core.w_states
import westpa.work_managers as work_managers
from common import MockArgs


class Test_W_States:
    def test_run_w_states(self, ref_initialized):
        '''Tests running w_states on a sample h5 file'''

        args = MockArgs(verbosity=0, rcfile=self.cfg_filepath)

        westpa.rc.process_args(args)
        work_managers.environment.process_wm_args(args)

        westpa.cli.core.w_states.initialize(mode='show', _bstate_file=None, bstates=None, tstates=None, _tstate_file=None)

        out, err = self.capfd.readouterr()

        state_string = out.split('\n')[-2]
        known_state_string = 'initial                1    None        # state_id=0    pcoord=[8.0]'

        assert state_string == known_state_string

    @fixture(autouse=True)
    def capfd(self, capfd):
        self.capfd = capfd
