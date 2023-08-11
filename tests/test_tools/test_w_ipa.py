import os
import shutil
from h5diff import H5Diff
import unittest
import pytest
import argparse
from unittest import mock

from westpa.cli.tools.w_ipa import entry_point


class Test_W_IPA(unittest.TestCase):
    test_name = 'W_IPA'

    def test_run_w_ipa(self):
        '''Testing if w_ipa runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_ref.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_ref.h5'), './west.h5')
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


@pytest.mark.skip(reason="work-in-progress test that uses entry point")
class Test_W_IPA_new:
    def test_run_w_ipa(self, ref_50iter):
        '''Testing if w_ipa runs as expected and the h5 files looks good.'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                rcfile=self.cfg_filepath,
                verbosity='debug',
                work_manager=None,
                analysis_mode=True,
                max_queue_length=None,
                debug_mode=True,
                we_h5filename=self.h5_filepath,
                scheme='TEST',
                reanalyze=False,
                ignore_hash=False,
                plotting=False,
                construct_dataset=False,
                dsspecs=None,
                output='assign.h5',
                subsample=None,
                config_from_file=True,
            ),
        ):
            entry_point()

        assert os.path.isfile('./ANALYSIS/TEST/assign.h5'), "The assign.h5 file was not generated."
        assert os.path.isfile('./ANALYSIS/TEST/direct.h5'), "The direct.h5 file was not generated."
        diff = H5Diff('assign_ipa_ref.h5', './ANALYSIS/TEST/assign.h5')
        diff.check()
        #       TODO: this is broken
        #       diff2 = H5Diff('direct_ipa_ref.h5', './ANALYSIS/TEST/direct.h5')
        #       diff2.check()

        # cleanup
        shutil.rmtree('ANALYSIS')
