import pytest
from filecmp import cmpfiles
from io import StringIO

import numpy as np
import pickle
from numpy.testing import assert_array_equal

import westpa
from westpa.core.propagators.executable import ExecutablePropagator
from westpa.core.propagators.loaders import (
    numpy_data_loader,
    pickle_data_loader,
    aux_data_loader,
    restart_loader,
    restart_writer,
    seglog_loader,
    seglog_writer,
    pcoord_loader,
)
from westpa.core.segment import Segment


class Test_Executable:
    """Class to test the propagator executable."""

    def test_data_config(self, ref_executable):
        """Test if the config is initialized correctly, where the executable propagator dataset options are set with the data manager options."""

        # Make the rc and executable read the config file.
        westpa.rc.read_config(filename='west_implicit.cfg')
        executable = ExecutablePropagator(rc=westpa.rc)

        assert 'displacement' in executable.data_info
        assert executable.data_info['displacement']['loader'] == numpy_data_loader

    def test_legacy_data_config(self, ref_executable):
        """Test if the dataset config is initialized correctly using the legacy part, where propagator datasets have to be specified twice."""

        # Make the rc and executable read the config file.
        westpa.rc.read_config(filename='west.cfg')
        executable = ExecutablePropagator(rc=westpa.rc)

        assert 'displacement' in executable.data_info
        assert executable.data_info['displacement']['loader'] == aux_data_loader


class Test_Loaders:
    """Class to test if numpy_data_loader and pickle_data_loader are able to successfully add data into a dummy segment object."""

    def test_numpy_loader(self, ref_idtype):
        """Test if data loaded with numpy_data_loader is consistent."""

        test_segment = Segment()

        numpy_data_loader('test', self.correct_pkl, test_segment, False)

        with open(self.correct_pkl, 'rb') as f:
            ref_array = pickle.load(f)

        test_array = test_segment.data['test'][:]

        assert np.array_equal(test_array, ref_array)

    def test_pickle_loader(self, ref_idtype):
        """Test if data loaded with numpy_data_loader is consistent."""

        test_segment = Segment()

        pickle_data_loader('test', self.correct_pkl, test_segment, False)

        with open(self.correct_pkl, 'rb') as f:
            ref_array = pickle.load(f)

        test_array = test_segment.data['test'][:]

        assert np.array_equal(test_array, ref_array)

    def test_restart_loader_writer(self, nacl_restart_files):
        """Test if the restart file can be read, saved and reloaded correctly."""

        # Make a dummy segment and read/write the restart files
        test_segment = Segment()
        restart_loader('restart', self.return_dir, test_segment, False)
        restart_writer(self.write_dir, test_segment)

        # Do a shallow file comparison and make sure files tarred up and written out matches
        (matches, mismatches, errors) = cmpfiles(self.return_dir, self.write_dir, ['nacl.prmtop', 'nacl.ncrst'])
        assert sum([True if file in self.nacl_restart_files else False for file in matches]) == 2

    def test_seglog_loader_writer(self, nacl_restart_files):
        """Test if the log file can be saved and reloaded correctly."""

        # Make a dummy segment and read/write the seglog file
        test_segment = Segment()

        # Generate dummy log file
        dummy_text = 'abc\nlog\nend'
        self.test_file_path = self.test_dir / 'test.log'
        with open(self.test_file_path, 'w') as text_file:
            text_file.write(dummy_text)

        # Save the file into
        seglog_loader('log', self.test_file_path, test_segment, False)

        # Write the current file into self.write_dir
        seglog_writer(self.write_dir, test_segment)

        # Check to ensure contents are preserved
        with open(self.write_dir / 'seg.log', 'r') as text_file:
            assert text_file.read() == dummy_text

    def test_pcoord_loader_failures(self, ref_mab):
        test_segment = Segment()

        # Making test data
        rng = np.random.default_rng()
        c = rng.random(size=(11, 2), dtype=np.float32)
        io_file = StringIO()
        np.savetxt(io_file, c)

        with pytest.raises(AssertionError):
            pcoord_loader('test', c, test_segment, False)

        with pytest.raises(ValueError, match=r'incorrect shape \(11, 2\) \[expected \(2, 1\)\]'):
            io_file.seek(0)
            pcoord_loader('pcoord', io_file, test_segment, False)

    def test_pcoord_loader(self, ref_mab):
        test_segment = Segment()

        # Making test data
        rng = np.random.default_rng()
        c = rng.random(size=(2, 1), dtype=np.float32)
        io_file = StringIO()
        np.savetxt(io_file, c)

        io_file.seek(0)
        pcoord_loader('pcoord', io_file, test_segment, False)

        assert_array_equal(test_segment.pcoord, c)
