import argparse

from .hdiff import H5Diff
import westpa

from westpa.cli.core.w_init import entry_point
from unittest import mock


class Test_W_Init:
    test_name = 'W_INIT'

    def test_run_w_init(self, initialized_west_ref):
        '''Tests initialization of a WESTPA simulation system from a prebuilt .cfg'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                force=True,
                rcfile=self.cfg_path,
                bstate_file=None,
                verbosity='0',
                bstates=['initial,1.0'],
                tstate_file=None,
                tstates=None,
                segs_per_state=1,
                shotgun=False,
            ),
        ):

            entry_point()

        # h5 files contain some internal information that includes timestamps, so I can't just compare md5 checksums
        #   to ensure that w_init is producing the same output.
        # Instead, use my H5Diff class.
        # If the checked contents differ, an AssertionError will be raised.
        diff = H5Diff(self.odld_path + '/west_ref.h5', self.odld_path + '/west.h5')
        diff.check()

    def tearDown(self):

        westpa.rc._sim_manager = None
        westpa.rc._system = None
        westpa.rc._data_manager = None
        westpa.rc._we_driver = None
        westpa.rc._propagator = None
