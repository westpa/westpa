import unittest
import argparse
import os

from .hdiff import H5Diff

from westpa.cli.core.w_run import entry_point
from unittest import mock

from shutil import copyfile


class Test_W_Run(unittest.TestCase):
    test_name = 'W_RUN'

    # This is a little kludgey, but in order to set the class attribute starting_path, I need to
    #   have an __init__ method or it errors out, despite the fact that I'm not actually initializing the
    #   variable in it
    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_run')

    def test_run_w_run(self):
        '''Tests running an initialized WESTPA system for 3 iterations'''

        self.starting_path = os.getcwd()

        odld_path = os.path.dirname(__file__) + '/ref'

        copyfile(odld_path + '/west_ref.h5', odld_path + '/west.h5')

        os.chdir(odld_path)

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(verbosity=0, rcfile=odld_path + '/west.cfg'),
        ):
            entry_point()

        # Verify that the generated h5 file is identical in contents to the 3-iteration reference
        diff = H5Diff(odld_path + '/west_3iter.h5', odld_path + '/west.h5')
        diff.check()

    def tearDown(self):

        os.remove(os.path.dirname(__file__) + '/odld/west.h5')
