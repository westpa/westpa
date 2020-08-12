import os
import shutil
import unittest
import filecmp


class Test_Plothist(unittest.TestCase):
    '''Class to test plothist function'''

    test_name = 'PLOTHIST'

    def test_compare_avg_txt(self):
        '''Test to compare plothist average reference text file to outputted text file'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'pdist_ref.h5'), './')
        shutil.copy2(os.path.join(ref_dir, 'hist_avg_ref.txt'), './')
        os.system("plothist average pdist_ref.h5 0 -o hist_avg.pdf --text-output phist_avg.txt --first-iter 1 --last-iter 50")
        os.system("cat phist_avg.txt | tail -n +10 > hist_avg.txt")
        assert os.path.isfile('./hist_avg.pdf'), "The average pdf was not generated."
        assert os.path.isfile('./hist_avg.txt'), "The text output file was not generated."
        self.assertTrue(filecmp.cmp('hist_avg_ref.txt', 'hist_avg.txt'), "Text output file not generated correctly.")
        os.remove('hist_avg.pdf')
        os.remove('phist_avg.txt')
        os.remove('hist_avg.txt')
        os.remove('hist_avg_ref.txt')
        os.remove('pdist_ref.h5')

    def test_compare_inst_txt(self):
        '''Test to compare plothist instant reference text file to outputted text file'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'pdist_ref.h5'), './')
        shutil.copy2(os.path.join(ref_dir, 'hist_inst_ref.txt'), './')
        os.system("plothist instant pdist_ref.h5 0 -o hist_inst.pdf --text-output phist_inst.txt --iter 50")
        os.system("cat phist_inst.txt | tail -n +10 > hist_inst.txt")
        assert os.path.isfile('./hist_inst.pdf'), "The instant pdf was not generated."
        assert os.path.isfile('./hist_inst.txt'), "The text output file was not generated."
        self.assertTrue(filecmp.cmp('hist_inst_ref.txt', 'hist_inst.txt'), "Text output file not generated correctly.")
        os.remove('hist_inst.pdf')
        os.remove('phist_inst.txt')
        os.remove('hist_inst.txt')
        os.remove('hist_inst_ref.txt')
        os.remove('pdist_ref.h5')

    def test_compare_evol_txt(self):
        '''Test to compare plothist evolution pdf was generated'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'pdist_ref.h5'), './')
        os.system("plothist evolution pdist_ref.h5 0 -o hist_evol.pdf --first-iter 1 --last-iter 50")
        assert os.path.isfile('./hist_evol.pdf'), "The evolution pdf was not generated."
        os.remove('hist_evol.pdf')
        os.remove('pdist_ref.h5')
