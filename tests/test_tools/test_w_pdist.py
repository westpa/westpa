import os
import shutil
from .h5diff import H5Diff
import unittest


class Test_W_PDIST(unittest.TestCase):
    '''Class to test w_pdist works to generate a file and that it is the same as the sample pdist.h5 file.'''

    test_name = 'W_PDIST'

    def test_run_w_pdist(self):
        '''Testing if w_pdist runs as expected and the pdist.h5 file looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'west.cfg'), './')
        shutil.copy2(os.path.join(ref_dir, 'west.h5'), './')
        os.system("w_pdist")
        assert os.path.isfile('./ANALYSIS/TEST/pdist.h5'), "The pdist.h5 file was not generated."
        diff = H5Diff(os.path.join(ref_dir, 'pdist_ref.h5'), './pdist.h5')
        diff.check()
        os.remove('west.h5')
        os.remove('pdist.h5')
