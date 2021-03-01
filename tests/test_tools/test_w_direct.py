import argparse
import os
import pytest
import shutil
import unittest

from h5diff import H5Diff

from westpa.cli.tools.w_direct import entry_point, DAll
from unittest import mock


class Test_W_Direct(unittest.TestCase):

    test_name = 'W_DIRECT'

    def test_run_w_direct(self):
        '''Testing if w_direct runs as expected and the direct.h5 file looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_ref.h5'), './west.h5')
        os.system('w_direct all')
        assert os.path.isfile('./direct.h5'), "The direct.h5 file was not generated."
        diff = H5Diff(os.path.join(ref_dir, 'direct_ref.h5'), './direct.h5')
        diff.check()
        os.remove('direct.h5')
        os.remove('west.h5')


@pytest.mark.skip(reason="work-in-progress test that uses entry point")
class Test_W_Direct_New:
    def test_run_w_direct(self, ref_50iter):

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile=self.cfg_filepath,
                work_manager=None,
                max_queue_length=None,
                west_subcommand=DAll(1),
                we_h5filename=self.h5_filepath,
                construct_dataset=None,
                dsspecs=None,
                output='direct.h5',
                kinetics='direct.h5',
                first_iter=1,
                last_iter=None,
                step_iter=None,
                assignments='assign_ref.h5',
                evolution_mode=None,
                subsample=None,
                config_from_file=True,
                scheme='TEST',
                bootstrap=None,
                correl=None,
                alpha=0.05,
                acalpha=None,
                nsets=None,
                window_frac=1.0,
                display_averages=True,
            ),
        ):
            entry_point()

        diff = H5Diff('./direct_ref.h5', './direct.h5')
        diff.check()

        # clean up
        os.remove('direct.h5')
