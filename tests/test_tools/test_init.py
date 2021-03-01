import argparse

from h5diff import H5Diff

import westpa.cli.core.w_init
from unittest import mock


class Test_W_Init:
    def test_run_w_init(self, ref_cfg):
        '''Tests initialization of a WESTPA simulation system from a prebuilt .cfg'''
        # This test is named in such a way so it always runs before test_w_assign.py. It will fail otherwise.

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                force=True,
                rcfile=self.cfg_filepath,
                bstate_file=None,
                verbosity='verbose',
                bstates=['initial,1.0'],
                tstate_file=None,
                tstates=None,
                segs_per_state=1,
                shotgun=False,
            ),
        ):
            westpa.cli.core.w_init.entry_point()

        # h5 files contain some internal information that includes timestamps, so I can't just compare md5 checksums
        #   to ensure that w_init is producing the same output.
        # Instead, use my H5Diff class.
        # If the checked contents differ, an AssertionError will be raised.
        diff = H5Diff(self.ref_h5_filepath, self.h5_filepath)
        diff.check()
