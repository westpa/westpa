import os
from unittest import mock
from h5diff import H5Diff
from westpa.cli.tools.w_multi_west import entry_point, create_idtype_array
import argparse
import numpy as np
import pickle
import h5py


class Test_W_Multi_West:
    '''Class to test w_pdist works to generate a file and that it is the same as the sample pdist.h5 file.'''

    def test_run_w_multi_west(self, ref_multi):
        '''Testing if w_multi_west runs as expected and the multi.h5 file looks good.'''

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
                ibstates=True,
            ),
        ):
            entry_point()

        assert os.path.isfile('./multi.h5'), "The multi.h5 file was not generated."
        diff = H5Diff('multi_aux_ref.h5', 'multi.h5')
        diff.check()

        # clean up
        os.remove('multi.h5')

    def test_run_w_multi_west_noaux(self, ref_multi_noaux):
        '''Testing if w_multi_west exception states runs as expected and the multi.h5 file looks good.'''

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
                ibstates=True,
            ),
        ):
            entry_point()

        assert os.path.isfile('./multi.h5'), "The multi.h5 file was not generated."
        diff = H5Diff('multi_noaux_ref.h5', 'multi.h5')
        diff.check()

        # clean up
        os.remove('multi.h5')

    def test_idtype(self, ref_idtype):
        '''Test if the new istate array is consistent if made with the create_idtype_array() function.'''

        input_array = h5py.File(self.h5_filepath)['ibstates/0/istate_index'][:]

        with open(self.correct_pkl, 'rb') as f:
            ref_array = pickle.load(f)

        test_array = create_idtype_array(input_array)

        assert np.array_equal(test_array, ref_array)
