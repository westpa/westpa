import os
import argparse
import numpy as np
from unittest import mock
from h5diff import H5Diff
from westpa.cli.tools.w_pdist import entry_point, WPDist


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
                construct_wdataset=None,
            ),
        ):
            entry_point()

        assert os.path.isfile('./pdist.h5'), "The pdist.h5 file was not generated."
        diff = H5Diff('pdist_ref.h5', 'pdist.h5')
        diff.check()

        # clean up
        os.remove('pdist.h5')

    def test_w_pdist_construct_bins(self):
        '''Test min/max autobinning'''
        correct_binbounds = np.asarray([[-20.0, -12.0, -4, 4.0, 12.0, 20.2] for _ in range(2)])
        correct_midpoint = np.asarray([range(-16, 17, 8) for _ in range(2)])

        pdist = WPDist()
        pdist.ndim = 2
        pdist.data_range = [(-20, 20), (-20, 20)]

        # Constructing with a single value
        pdist.construct_bins(5)
        print(f'{pdist.midpoints=}')
        assert np.array_equal(pdist.binbounds, correct_binbounds)
        assert np.array_equal(pdist.midpoints, correct_midpoint)

        # Constructing with a list of bin counts
        pdist.construct_bins([5, 5])
        assert np.array_equal(pdist.binbounds, correct_binbounds)
        assert np.array_equal(pdist.midpoints, correct_midpoint)

        pdist.construct_bins(np.asarray([range(-20, 21, 8) for _ in range(2)]))
        print(f'{pdist.binbounds=}')
        print(f'{correct_binbounds=}')
        assert np.array_equal(pdist.binbounds, correct_binbounds)
        assert np.array_equal(pdist.midpoints, correct_midpoint)
