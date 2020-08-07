import sys, os, shutil, tempfile

from h5diff import H5Diff

import unittest

import h5py

sys.dont_write_bytecode = True

class Test_W_IPA(unittest.TestCase):

    test_name = 'W_IPA'

    def test_run_w_ipa(self):
        '''Testing if w_ipa runs as expected and the h5 files looks good.'''

        temp = tempfile.TemporaryDirectory(dir = "./")
        os.chdir(temp.name)
        shutil.copy2('../refs/west.cfg', './')
        shutil.copy2('../refs/west.h5', './')
        os.system("w_ipa -ao")
        assert os.path.isfile('./ANALYSIS/TEST/assign.h5'), "The assign.h5 file was not generated."
        assert os.path.isfile('./ANALYSIS/TEST/direct.h5'), "The direct.h5 file was not generated."

        diff = H5Diff('../refs/assign_ipa_ref.h5', './ANALYSIS/TEST/assign.h5')
#       TODO: this is broken
#       diff = H5Diff('../refs/direct_ipa_ref.h5', './ANALYSIS/TEST/direct.h5')
        diff.check()

        os.chdir("../")
        temp.cleanup()
