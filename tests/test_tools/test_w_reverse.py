import argparse
import os
import shutil
import unittest
from filecmp import cmp
from westpa.cli.tools.w_reverse import entry_point


class Test_W_Reverse(unittest.TestCase):
    test_name = 'W_Reverse'

    def test_run_w_reverse_no_hdf5(self):
        '''Testing if w_reverse runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_no_hdf5.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_no_hdf5.h5'), './west.h5')
        shutil.copytree(os.path.join(ref_dir, 'traj_segs_reverse_no_hdf5'), './traj_segs')
        with unittest.mock.patch(
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
                rcfile='west.cfg',
                verbosity='debug',
            ),
        ):
            entry_point()
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert os.path.isfile('./bstates_reverse/000001_000000.xml'), "The 000001_000000.xml file was not generated."
        assert os.path.getsize('./bstates_reverse/000001_000000.xml') > 0, "The 000001_000000.xml file is empty."
        assert os.path.isfile('./bstates_reverse/000002_000000.xml'), "The 000002_000000.xml file was not generated."
        assert os.path.getsize('./bstates_reverse/000002_000000.xml') > 0, "The 000002_000000.xml file is empty."
        assert os.path.isfile('./bstates_reverse/000003_000000.xml'), "The 000003_000000.xml file was not generated."
        assert os.path.getsize('./bstates_reverse/000003_000000.xml') > 0, "The 000003_000000.xml file is empty."
        assert cmp(
            os.path.join(ref_dir, 'bstates.txt'), './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are not the same'
        shutil.rmtree('traj_segs')
        os.remove('west.h5')
        os.remove('west.cfg')

    def test_run_w_reverse_hdf5_no_rst_file(self):
        '''Testing if w_reverse runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.h5'), './west.h5')
        shutil.copytree(os.path.join(ref_dir, 'traj_segs_reverse_hdf5'), './traj_segs')
        with unittest.mock.patch(
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
                rcfile='west.cfg',
                verbosity='debug',
            ),
        ):
            entry_point()
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        assert os.path.isfile('./bstates_reverse/000001_000000.xml'), "The 000001_000000.xml file was not generated."
        assert os.path.getsize('./bstates_reverse/000001_000000.xml') > 0, "The 000001_000000.xml file is empty."
        assert os.path.isfile('./bstates_reverse/000002_000000.xml'), "The 000002_000000.xml file was not generated."
        assert os.path.getsize('./bstates_reverse/000002_000000.xml') > 0, "The 000002_000000.xml file is empty."
        assert os.path.isfile('./bstates_reverse/000003_000000.xml'), "The 000003_000000.xml file was not generated."
        assert os.path.getsize('./bstates_reverse/000003_000000.xml') > 0, "The 000003_000000.xml file is empty."
        assert cmp(
            os.path.join(ref_dir, 'bstates.txt'), './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are not the same'
        shutil.rmtree('traj_segs')
        os.remove('west.h5')
        os.remove('west.cfg')

    def test_run_w_reverse_hdf5_rst_file(self):
        '''Testing if w_reverse runs as expected and the h5 files looks good.'''

        ref_dir = os.path.join(os.path.dirname(__file__), '../refs')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.cfg'), './west.cfg')
        shutil.copy2(os.path.join(ref_dir, 'west_reverse_hdf5.h5'), './west.h5')
        shutil.copytree(os.path.join(ref_dir, 'traj_segs_reverse_hdf5'), './traj_segs')
        with unittest.mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                we_h5filename='west.h5',
                first_iter=1,
                last_iter=None,
                config_file='west.cfg',
                max_n_bstates=10000,
                rst_file='parent.xml',
                output_bstates_dir='bstates_reverse',
                output_bstates_file='bstates.txt',
                use_weights=True,
                seed=12345,
                rcfile='west.cfg',
                verbosity='debug',
            ),
        ):
            entry_point()
        assert os.path.isfile('./bstates_reverse/bstates.txt'), "The bstates.txt file was not generated."
        for iiter in range(1, 4):
            assert os.path.isfile(f'./bstates_reverse/{iiter:06d}_000000.xml'), "The {iiter:06d}_000000.xml file was not generated."
            assert os.path.getsize(f'./bstates_reverse/{iiter:06d}_000000.xml') > 0, f"The {iiter:06d}_000000.xml file is empty."
        assert cmp(
            os.path.join(ref_dir, 'bstates.txt'), './bstates_reverse/bstates.txt'
        ), 'The reference bstates.txt and the produced bstates.txt are not the same'
        shutil.rmtree('traj_segs')
        os.remove('west.h5')
        os.remove('west.cfg')
