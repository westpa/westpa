import os
import shutil
from .h5diff import H5Diff
import unittest


class Test_W_Assign(unittest.TestCase):

    test_name = 'W_ASSIGN'

    def test_run_w_assign(self):
        '''Testing if w_assign runs as expected and the assign.h5 file looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'west.cfg'), './')
        shutil.copy2(os.path.join(ref_dir, 'west.h5'), './')
        os.system('w_assign -W ./west.h5 --config-from-file --scheme TEST')
        assert os.path.isfile('./ANALYSIS/TEST/assign.h5'), "The assign.h5 file was not generated."
        diff = H5Diff(os.path.join(ref_dir, 'assign_ref.h5'), './ANALYSIS/TEST/assign.h5')
        diff.check()

        shutil.rmtree('ANALYSIS')
        os.remove('west.h5')
        os.remove('west.cfg')
