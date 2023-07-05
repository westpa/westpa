from .test_tools.conftest import *  # noqa
from westpa.core.propagators.executable import npy_data_loader, pickle_data_loader
from westpa.core.segment import Segment

import numpy as np
import pickle


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
