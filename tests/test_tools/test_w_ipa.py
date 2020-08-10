import os
import shutil
from .h5diff import H5Diff
import unittest


class Test_W_IPA(unittest.TestCase):

    test_name = 'W_IPA'

    def test_run_w_ipa(self):
        '''Testing if w_ipa runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'west.cfg'), './')
        shutil.copy2(os.path.join(ref_dir, 'west.h5'), './')
        os.system('w_ipa -ao')
        assert os.path.isfile('./ANALYSIS/TEST/assign.h5'), "The assign.h5 file was not generated."
        assert os.path.isfile('./ANALYSIS/TEST/direct.h5'), "The direct.h5 file was not generated."
        diff = H5Diff(os.path.join(ref_dir, 'assign_ipa_ref.h5'), './ANALYSIS/TEST/assign.h5')
        #       TODO: this is broken
        #       diff = H5Diff('../refs/direct_ipa_ref.h5', './ANALYSIS/TEST/direct.h5')
        diff.check()
        shutil.rmtree('ANALYSIS')
        os.remove('west.h5')
        os.remove('west.cfg')
