import h5py

import westpa.cli.core.w_run
import westpa.work_managers as work_managers
from common import MockArgs


class Test_W_Run:
    def test_run_w_run(self, ref_initialized):
        '''Tests running an initialized WESTPA system for 3 iterations'''

        args = MockArgs(
            rcfile=self.cfg_filepath,
            verbosity='0',
        )

        westpa.rc.process_args(args)
        work_managers.environment.process_wm_args(args)

        westpa.cli.core.w_run.run_simulation()

        # Since the dynamics produce slightly different output in each iteration, just check that
        #    the third iteration exists in the h5 file.
        # TODO: Modify the ODLD system to run identically across runs, and then remove the excluded
        #   datasets.

        _hfile = h5py.File(self.h5_filepath)
        assert 'iter_00000003' in _hfile['/iterations'].keys()
        _hfile.close()
