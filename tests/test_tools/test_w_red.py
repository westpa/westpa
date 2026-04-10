from common import MockArgs

import h5py
import westpa
from westpa.cli.tools import w_red


class Test_W_RED:
    '''Class to test w_red works to generate a dataset in direct.h5 file.'''

    def test_run_w_red(self, ref_50iter):
        '''Testing if w_red runs as expected and the `red_flux_evolution` dataset is made within direct.h5.'''

        rc = westpa.rc

        args = MockArgs(
            verbosity='debug',
            rcfile=self.cfg_filepath,
            max_queue_length=None,
            we_h5filename=self.h5_filepath,
            compress=False,
            work_manager=None,
            n_workers=None,
        )

        rc.process_args(args)

        tool = w_red.WRed()

        # Prepare and instantiate work manager
        tool.wm_env.process_wm_args(args)
        tool.work_manager = tool.wm_env.make_work_manager()

        tool.process_all_args(args)
        with tool.work_manager:
            if tool.work_manager.is_master:
                tool.go()
            else:
                tool.work_manager.run()

        with h5py.File('ANALYSIS/TEST/direct.h5', 'r') as h5file:
            assert 'red_flux_evolution' in h5file
