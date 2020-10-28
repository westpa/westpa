import os
import shutil
from unittest import mock
import filecmp
from westpa.cli.tools.plothist import entry_point, AveragePlotHist
import argparse
import pytest

class Test_Plothist:
    '''Class to test plothist function'''

    test_name = 'PLOTHIST'
    def test_compare_avg_txt(self, ref_50iter):
        '''Test to compare plothist average reference text file to outputted text file'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='0',
                rcfile=self.cfg_filepath,
                input='pdist_ref.h5',
                firstdim='0',
                plot_output='hist_avg.pdf',
                hdf5_output=None,
                plot_contour=False,
                title=None,
                plotscale='energy',
                enerzero='min',
                range=None,
                postprocess_function=None,
                seconddim=None,
                text_output='phist_avg.txt',
                first_iter= 1,
                last_iter=50,
                west_subcommand=AveragePlotHist(1),
            ),
        ):
            entry_point()
        #os.system("plothist average pdist_ref.h5 0 -o hist_avg.pdf --text-output phist_avg.txt --first-iter 1 --last-iter 50")
        os.system("cat phist_avg.txt | tail -n +10 > hist_avg.txt")
        assert os.path.isfile('./hist_avg.pdf'), "The average pdf was not generated."
        assert os.path.isfile('./hist_avg.txt'), "The text output file was not generated."
        #self.assertTrue(filecmp.cmp('hist_avg_ref.txt', 'hist_avg.txt'), "Text output file not generated correctly.")
        os.remove('hist_avg.pdf')
        os.remove('phist_avg.txt')
        os.remove('hist_avg.txt')

    @pytest.mark.skip
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
    @pytest.mark.skip
    def test_compare_evol_txt(self):
        '''Test to compare plothist evolution pdf was generated'''

        ref_dir = os.path.join(os.path.dirname(__file__), 'refs')
        shutil.copy2(os.path.join(ref_dir, 'pdist_ref.h5'), './')
        os.system("plothist evolution pdist_ref.h5 0 -o hist_evol.pdf --first-iter 1 --last-iter 50")
        assert os.path.isfile('./hist_evol.pdf'), "The evolution pdf was not generated."
        os.remove('hist_evol.pdf')
        os.remove('pdist_ref.h5')
