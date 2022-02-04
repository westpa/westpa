import os
from unittest import mock
from h5diff import H5Diff
from westpa.cli.tools.w_multi_west import entry_point
import argparse


class Test_W_Multi_West:
    '''Class to test w_pdist works to generate a file and that it is the same as the sample pdist.h5 file.'''

    def test_run_w_multi_west(self, ref_multi):
        '''Testing if w_pdist runs as expected and the pdist.h5 file looks good.'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                verbosity='debug',
                rcfile=self.cfg_filepath,
                max_queue_length=None,
                west=self.h5_filepath,
                output_file='multi.h5',
                master='.',
                sims='3',
                aux=None,
                auxall=True,
                no_reweight=False,
            ),
        ):
            entry_point()

        assert os.path.isfile('./multi.h5'), "The multi.h5 file was not generated."
        diff = H5Diff('multi_aux_ref.h5', 'multi.h5')
        diff.check()

        # clean up
        os.remove('multi.h5')
