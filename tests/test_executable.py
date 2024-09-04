import numpy as np
import pickle
import westpa

from westpa.core.propagators.executable import npy_data_loader, pickle_data_loader, aux_data_loader, ExecutablePropagator
from westpa.core.segment import Segment


class Test_Executable:
    '''Class to test the propagator executable.'''

    def test_data_config(self, ref_executable):
        '''Test if the config is initialized correctly, where the executable propagator dataset options are set with the data manager options.'''

        # Make the rc and executable read the config file.
        westpa.rc.read_config(filename='west_implicit.cfg')
        executable = ExecutablePropagator(rc=westpa.rc)

        assert 'displacement' in executable.data_info
        assert executable.data_info['displacement']['loader'] == npy_data_loader

    def test_legacy_data_config(self, ref_executable):
        '''Test if the dataset config is initialized correctly using the legacy part, where propagator datasets have to be specified twice.'''

        # Make the rc and executable read the config file.
        westpa.rc.read_config(filename='west.cfg')
        executable = ExecutablePropagator(rc=westpa.rc)

        assert 'displacement' in executable.data_info
        assert executable.data_info['displacement']['loader'] == aux_data_loader


class Test_Loaders:
    '''Class to test if npy_data_loader and pickle_date_loader are able to successfully add data into a dummy segment object.'''

    def test_npy_loader(self, ref_idtype):
        '''Test if data loaded with npy_data_loader is consistent.'''

        test_segment = Segment()

        npy_data_loader('test', self.correct_pkl, test_segment, False)

        with open(self.correct_pkl, 'rb') as f:
            ref_array = pickle.load(f)

        test_array = test_segment.data['test'][:]

        assert np.array_equal(test_array, ref_array)

    def test_pickle_loader(self, ref_idtype):
        '''Test if data loaded with npy_data_loader is consistent.'''

        test_segment = Segment()

        pickle_data_loader('test', self.correct_pkl, test_segment, False)

        with open(self.correct_pkl, 'rb') as f:
            ref_array = pickle.load(f)

        test_array = test_segment.data['test'][:]

        assert np.array_equal(test_array, ref_array)
