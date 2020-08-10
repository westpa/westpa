import os
import shutil
from .h5diff import H5Diff
import unittest


class Test_W_Assign(unittest.TestCase):

    test_name = 'W_DIRECT'

    def test_run_w_direct(self):
        '''Testing if w_direct runs as expected and the direct.h5 file looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'west.h5'), './')
        os.system('w_direct all')
        assert os.path.isfile('./direct.h5'), "The direct.h5 file was not generated."
        diff = H5Diff(os.path.join(ref_dir, 'direct_ref.h5'), './direct.h5')
        diff.check()
        os.remove('direct.h5')
        os.remove('west.h5')
