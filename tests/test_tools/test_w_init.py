import unittest
import argparse
import os

from .common import CommonToolTest
from .hdiff import H5Diff

from westpa.cli.core.w_init import entry_point
from unittest import mock


class Test_W_Init(unittest.TestCase, CommonToolTest):
    def test_init(self):
        '''Tests initialization of a WESTPA simulation system from a prebuilt .cfg'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                force=True,
                rcfile='west.cfg',
                bstate_file=None,
                verbosity=None,
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
        diff = H5Diff('west_ref.h5', 'west.h5')
        diff.check()

    def tearDown(self):

        # Don't need this any more
        os.remove('west.h5')
