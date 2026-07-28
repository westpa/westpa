import os
import shutil
import unittest
from filecmp import cmp


class Test_W_Reverse(unittest.TestCase):
    test_name = 'W_Reverse'

    def test_run_w_reverse_no_hdf5(self):
        '''Testing if w_reverse runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_no_hdf5.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_no_hdf5.h5'), './west.h5')
        shutil.copytree(os.path.join(ref_dir, 'traj_segs_reverse_no_hdf5'), './traj_segs')
        os.system('w_reverse --rst-file seg.xml')
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert os.path.isfile('./bstates_reverse/000001_000000.xml'), "The 000001_000000.xml file was not generated."
        assert os.path.isfile('./bstates_reverse/000002_000000.xml'), "The 000002_000000.xml file was not generated."
        assert os.path.isfile('./bstates_reverse/000003_000000.xml'), "The 000003_000000.xml file was not generated."
        assert cmp(
            os.path.join(ref_dir, 'bstates_no_hdf5.txt'), './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are notthe same'
        shutil.rmtree('traj_segs')
        os.remove('west.h5')
        os.remove('west.cfg')

    def test_run_w_reverse_hdf5_no_rst_file(self):
        '''Testing if w_reverse runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.h5'), './west.h5')
        shutil.copytree(os.path.join(ref_dir, 'traj_segs_reverse_hdf5'), './traj_segs')
        os.system('w_reverse')
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert os.path.isfile('./bstates_reverse/000002_000000.xml'), "The 000002_000000.xml file was not generated."
        assert os.path.isfile('./bstates_reverse/000003_000000.xml'), "The 000003_000000.xml file was not generated."
        assert cmp(
            os.path.join(ref_dir, 'bstates_hdf5.txt'), './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are notthe same'
        shutil.rmtree('traj_segs')
        os.remove('west.h5')
        os.remove('west.cfg')

    def test_run_w_reverse_hdf5_rst_file(self):
        '''Testing if w_reverse runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.h5'), './west.h5')
        shutil.copytree(os.path.join(ref_dir, 'traj_segs_reverse_hdf5'), './traj_segs')
        os.system('w_reverse --rst-file parent.xml')
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert os.path.isfile('./bstates_reverse/000002_000000.xml'), "The 000002_000000.xml file was not generated."
        assert os.path.isfile('./bstates_reverse/000003_000000.xml'), "The 000003_000000.xml file was not generated."
        assert cmp(
            os.path.join(ref_dir, 'bstates_hdf5.txt'), './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are notthe same'
        shutil.rmtree('traj_segs')
        os.remove('west.h5')
        os.remove('west.cfg')
