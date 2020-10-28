import argparse
import os
import pytest

from h5diff import H5Diff

from westpa.cli.tools.w_direct import entry_point, DAll
from unittest import mock


@pytest.mark.skip
class Test_W_Direct:
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
