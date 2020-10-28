import argparse
import shutil

from h5diff import H5Diff

from westpa.cli.tools.w_assign import entry_point
from unittest import mock


class Test_W_Assign:
    def test_run_w_assign(self, ref_50iter):

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile=self.cfg_filepath,
                max_queue_length=None,
                we_h5filename=self.h5_filepath,
                construct_dataset=None,
                dsspecs=None,
                output='assign.h5',
                subsample=None,
                config_from_file=True,
                scheme='TEST',
            ),
        ):
            entry_point()

        diff = H5Diff('./assign_ref.h5', './ANALYSIS/TEST/assign.h5')
        diff.check()

        # clean up
        shutil.rmtree('ANALYSIS')
