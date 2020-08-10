import unittest
import argparse
import os

from .hdiff import H5Diff
import westpa

from westpa.cli.core.w_init import entry_point
from unittest import mock


class Test_W_Init(unittest.TestCase):
    test_name = 'W_INIT'

    # This is a little kludgey, but in order to set the class attribute starting_path in runTest, I need to
    #   have an __init__ method or it errors out, despite the fact that I don't actually initialize the
    #   variable in it
    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_init')

    def setUp(self):

        # TODO: This *shouldn't* need to be tracked, after calling chdir everything should be dumped into odld_path
        #    However, it seems that this is not the case, and the generated west.h5 is dropped wherever pytest is
        #    launched from, hence the need for tracking this.
        self.starting_path = os.getcwd()

        self.odld_path = os.path.dirname(__file__) + '/ref'
        os.environ['WEST_SIM_ROOT'] = self.odld_path

        os.chdir(self.odld_path)

    def test_run_w_init(self):
        '''Tests initialization of a WESTPA simulation system from a prebuilt .cfg'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                force=True,
                rcfile=self.odld_path + '/west.cfg',
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
        os.environ['WEST_SIM_ROOT'] = ''
        os.remove(self.odld_path + '/west.h5')
        os.chdir(self.starting_path)
