import os
import shutil
import unittest


class Test_W_Reverse(unittest.TestCase):
    test_name = 'W_Reverse'

    def test_run_w_reverse(self):
        '''Testing if w_reverse runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse.h5'), './west.h5')
        shutil.copytree(os.path.join(ref_dir, 'traj_segs_reverse'), './traj_segs')
        os.system('w_reverse --rst-file seg.xml')
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert os.path.isfile('./bstates_reverse/000001_000000.xml'), "The 000001_000000.xml file was not generated."
        assert os.path.isfile('./bstates_reverse/000002_000000.xml'), "The 000002_000000.xml file was not generated."
        assert os.path.isfile('./bstates_reverse/000003_000000.xml'), "The 000003_000000.xml file was not generated."
        same_lines = True
        with open(os.path.join(ref_dir, 'bstates.txt'), 'r') as ref_file, open('./bstates_reverse/bstates.txt', 'r') as test_file:
            ref_lines = ref_file.readlines()
            test_lines = test_file.readlines()
            for line in test_lines:
                if line not in ref_lines:
                    same_lines = False
                    break
        assert same_lines, 'The reference bstates.txt and the produced bstates.txt do not contain the same information'
        shutil.rmtree('traj_segs')
        os.remove('west.h5')
        os.remove('west.cfg')
