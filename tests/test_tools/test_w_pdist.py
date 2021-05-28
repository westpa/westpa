import os
from unittest import mock
from h5diff import H5Diff
from westpa.cli.tools.w_pdist import entry_point
import argparse


class Test_W_PDIST:
    '''Class to test w_pdist works to generate a file and that it is the same as the sample pdist.h5 file.'''

    def test_run_w_pdist(self, ref_50iter):
        '''Testing if w_pdist runs as expected and the pdist.h5 file looks good.'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile=self.cfg_filepath,
                max_queue_length=None,
                we_h5filename=self.h5_filepath,
                construct_dataset=None,
                dsspecs=None,
                group_name='pcoord',
                first_iter=1,
                last_iter=None,
                bins='100',
                output='pdist.h5',
                ignore_out_of_range=False,
                compress=False,
                work_manager=None,
                n_workers=None,
            ),
        ):
            entry_point()

        assert os.path.isfile('./pdist.h5'), "The pdist.h5 file was not generated."
        diff = H5Diff('pdist_ref.h5', 'pdist.h5')
        diff.check()

        # clean up
        os.remove('pdist.h5')
