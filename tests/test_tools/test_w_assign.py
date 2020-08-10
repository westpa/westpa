import sys, os, shutil, tempfile

from .h5diff import H5Diff

import unittest

import h5py

sys.dont_write_bytecode = True

class Test_W_Assign(unittest.TestCase):

    test_name = 'W_ASSIGN'

    def test_run_w_assign(self):
        '''Testing if w_assign runs as expected and the assign.h5 file looks good.'''

        temp = tempfile.TemporaryDirectory(dir = "./")
        os.chdir(temp.name)
        shutil.copy2('../refs/west.cfg', './')
        shutil.copy2('../refs/west.h5', './')
        os.system("w_assign --config-from-file --scheme TEST")
        assert os.path.isfile('./ANALYSIS/TEST/assign.h5'), "The assign.h5 file was not generated."

        diff = H5Diff('../refs/assign_ref.h5', './ANALYSIS/TEST/assign.h5')
        diff.check()

        os.chdir("../")
        temp.cleanup()
