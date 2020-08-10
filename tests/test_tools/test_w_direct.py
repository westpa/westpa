import sys, os, shutil, tempfile

from .h5diff import H5Diff

import unittest

import h5py

sys.dont_write_bytecode = True

class Test_W_Assign(unittest.TestCase):

    test_name = 'W_DIRECT'

    def test_run_w_direct(self):
        '''Testing if w_direct runs as expected and the direct.h5 file looks good.'''

        temp = tempfile.TemporaryDirectory(dir = "./")
        os.chdir(temp.name)
        shutil.copy2('../refs/west.h5', './')
        os.system("w_direct all")
        assert os.path.isfile('./direct.h5'), "The direct.h5 file was not generated."

        diff = H5Diff('../refs/direct_ref.h5', './direct.h5')
        diff.check()

        os.chdir("../")
        temp.cleanup()
