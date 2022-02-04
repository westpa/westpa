from h5diff import H5Diff
import westpa.cli.core.w_init
from common import MockArgs


class Test_W_Init:
    def test_run_w_init(self, ref_cfg):
        '''Tests initialization of a WESTPA simulation system from a prebuilt .cfg'''
        # This test is named in such a way so it always runs before test_w_assign.py. It will fail otherwise.

        # The argument processing is just looking for an object with these attributes
        args = MockArgs(force=True, rcfile=self.cfg_filepath, verbosity='verbose')

        westpa.rc.process_args(args)

        westpa.cli.core.w_init.initialize(
            tstates=None,
            tstate_file=None,
            bstates=['initial,1.0'],
            bstate_file=None,
            sstates=None,
            sstate_file=None,
            segs_per_state=1,
            shotgun=False,
        )

        # h5 files contain some internal information that includes timestamps, so I can't just compare md5 checksums
        #   to ensure that w_init is producing the same output.
        # Instead, use my H5Diff class.
        # If the checked contents differ, an AssertionError will be raised.
        diff = H5Diff(self.ref_h5_filepath, self.h5_filepath)
        diff.check()
