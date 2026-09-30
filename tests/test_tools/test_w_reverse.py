import argparse
import os
from filecmp import cmp
from unittest import mock

import numpy as np

from westpa.cli.tools.w_reverse import W_Reverse, entry_point


class Test_W_Reverse:
    test_name = 'W_Reverse'

    def test_run_w_reverse_no_hdf5(self, w_reverse_bstate_no_hdf5_files):
        '''Testing if w_reverse (w/o hdf5 framework) runs as expected and resulting files looks good.'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                we_h5filename='west.h5',
                first_iter=1,
                last_iter=None,
                config_file='west.cfg',
                max_n_bstates=10000,
                rst_file='seg.xml',
                output_bstates_dir='bstates_reverse',
                output_bstates_file='bstates.txt',
                use_weights=True,
                seed=12345,
                verbosity=10,
                rcfile='west.cfg',
            ),
        ):
            entry_point()

        for iiter in range(1, 4):
            assert os.path.isfile(f'./bstates_reverse/{iiter:06d}_000000.xml'), "The {iiter:06d}_000000.xml file was not generated."
            assert os.path.getsize(f'./bstates_reverse/{iiter:06d}_000000.xml'), f"The {iiter:06d}_000000.xml file is empty."
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert cmp(
            'bstates.txt', './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are not the same'

    def test_run_w_reverse_hdf5_no_rst_file(self, w_reverse_bstate_hdf5_files):
        '''Testing if w_reverse (w/ hdf5 framework) runs as expected with `rst_file=None` and the resulting files looks good.'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                we_h5filename='west.h5',
                first_iter=1,
                last_iter=None,
                config_file='west.cfg',
                max_n_bstates=10000,
                rst_file=None,
                output_bstates_dir='bstates_reverse',
                output_bstates_file='bstates.txt',
                use_weights=True,
                seed=12345,
                verbosity=10,
                rcfile='west.cfg',
            ),
        ):
            entry_point()

        for iiter in range(1, 4):
            assert os.path.isfile(f'./bstates_reverse/{iiter:06d}_000000.xml'), "The {iiter:06d}_000000.xml file was not generated."
            assert os.path.getsize(f'./bstates_reverse/{iiter:06d}_000000.xml'), f"The {iiter:06d}_000000.xml file is empty."
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert cmp(
            'bstates.txt', './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are not the same'

    def test_run_w_reverse_hdf5_rst_file(self, w_reverse_bstate_hdf5_files):
        '''Testing if w_reverse (w/ hdf5 framework) runs as expected with `rst_file` defined `and the h5 files looks good.'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                we_h5filename='west.h5',
                first_iter=1,
                last_iter=None,
                config_file='west.cfg',
                max_n_bstates=10000,
                rst_file='seg.xml',
                output_bstates_dir='bstates_reverse',
                output_bstates_file='bstates.txt',
                use_weights=True,
                seed=12345,
                verbosity=10,
                rcfile='west.cfg',
            ),
        ):
            entry_point()

        for iiter in range(1, 4):
            assert os.path.isfile(f'./bstates_reverse/{iiter:06d}_000000.xml'), "The {iiter:06d}_000000.xml file was not generated."
            assert os.path.getsize(f'./bstates_reverse/{iiter:06d}_000000.xml'), f"The {iiter:06d}_000000.xml file is empty."
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert cmp(
            'bstates.txt', './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are not the same'

    def test_find_recycled_segments(self, w_reverse_bstate_hdf5_files):
        wr = W_Reverse()
        wr.data_reader.we_h5filename = 'west.h5'
        with wr.data_reader:
            succ_pairs = wr._find_recycled_segments()

        np.array_equal(succ_pairs, [(0, 0, 1.0), (1, 0, 1.0), (2, 0, 1.0)])
