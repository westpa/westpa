import os
import shutil
import tempfile
import filecmp
import unittest

# from westpa.cli.core.w_succ import find_successful_trajs


class Test_W_Succ(unittest.TestCase):
    '''
    Class to test that w_succ is working successfully to find and output recycled trajectories.
    '''

    test_name = 'W_SUCC'

    def test_w_succ_output(self):
        '''
        First, test if the output contains the correct information of segments which successfully reach a target state.
        Then test for correct file output with -o and --output flags.
        '''

        ref_dir = os.path.join(os.getcwd(), "refs")
        w_succ_ref = os.path.join(ref_dir, "w_succ_ref")
        temp = tempfile.TemporaryDirectory(dir=os.getcwd())

        os.chdir(temp.name)
        shutil.copy2(os.path.join(ref_dir, "west.cfg"), temp.name)
        shutil.copy2(os.path.join(ref_dir, "west.h5"), temp.name)

        # generate w_succ command line output file and standardize it by removing WEST data dir header line
        os.system("w_succ > temp.txt")
        os.system("sed -i'.original' '/#/,$!d' temp.txt")

        self.assertTrue(filecmp.cmp("temp.txt", w_succ_ref), "the cli output was not successfully generated")

        # test output to file args
        os.system("w_succ --output temp_test_output.txt")
        self.assertTrue(filecmp.cmp("temp_test_output.txt", w_succ_ref), "--output flag output file was not successfully generated")

        os.system("w_succ -o temp_test_o.txt")
        self.assertTrue(filecmp.cmp("temp_test_o.txt", w_succ_ref), "-o flag output file was not successfully generated")

        os.chdir(os.path.dirname(os.getcwd()))
        temp.cleanup()
