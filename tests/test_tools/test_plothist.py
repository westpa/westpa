import argparse
import os
from unittest import mock

import numpy as np

from h5diff import H5Diff
from westpa.cli.tools.plothist import entry_point, AveragePlotHist, InstantPlotHist, EvolutionPlotHist


class Test_Plothist:
    '''Class to test plothist function'''

    test_name = 'PLOTHIST'

    def test_compare_avg_txt(self, ref_50iter):
        '''Test to compare plothist average reference text file to outputted text file'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
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
                first_iter=1,
                last_iter=50,
                west_subcommand=AveragePlotHist(1),
            ),
        ):
            entry_point()

        os.system("cat phist_avg.txt | tail -n +10 > hist_avg.txt")
        assert os.path.isfile('./hist_avg.pdf'), "The average pdf was not generated."
        assert os.path.isfile('./hist_avg.txt'), "The text output file was not generated."

        expected = np.loadtxt('hist_avg_ref.txt')
        actual = np.loadtxt('hist_avg.txt')
        assert np.allclose(actual, expected), "Text output file not generated correctly."

        os.remove('hist_avg.pdf')
        os.remove('phist_avg.txt')
        os.remove('hist_avg.txt')

    def test_compare_inst_txt(self, ref_50iter):
        '''Test to compare plothist instant reference text file to outputted text file'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile=self.cfg_filepath,
                input='pdist_ref.h5',
                firstdim='0',
                plot_output='hist_inst.pdf',
                hdf5_output=None,
                plot_contour=False,
                title=None,
                plotscale='energy',
                enerzero='min',
                range=None,
                postprocess_function=None,
                seconddim=None,
                text_output='phist_inst.txt',
                first_iter=1,
                n_iter=50,
                west_subcommand=InstantPlotHist(1),
            ),
        ):
            entry_point()

        os.system("cat phist_inst.txt | tail -n +10 > hist_inst.txt")
        assert os.path.isfile('./hist_inst.pdf'), "The instant pdf was not generated."
        assert os.path.isfile('./hist_inst.txt'), "The text output file was not generated."

        expected = np.loadtxt('hist_inst_ref.txt')
        actual = np.loadtxt('hist_inst.txt')
        assert np.allclose(actual, expected), "Text output file not generated correctly."

        os.remove('hist_inst.pdf')
        os.remove('phist_inst.txt')
        os.remove('hist_inst.txt')

    def test_compare_evol_h5(self, ref_50iter):
        '''Test to compare plothist evolution pdf was generated'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile=self.cfg_filepath,
                input='pdist_ref.h5',
                firstdim='0',
                plot_output='hist_evol.pdf',
                hdf5_output='hist_evol.h5',
                plot_contour=False,
                title=None,
                plotscale='energy',
                enerzero='min',
                range=None,
                postprocess_function=None,
                seconddim=None,
                first_iter=1,
                last_iter=50,
                step_iter=None,
                west_subcommand=EvolutionPlotHist(1),
            ),
        ):
            entry_point()

        assert os.path.isfile('./hist_evol.pdf'), "The evolution pdf was not generated."
        diff = H5Diff('hist_evol_ref.h5', 'hist_evol.h5')
        diff.check()

        # clean up
        os.remove('hist_evol.pdf')
        os.remove('hist_evol.h5')
