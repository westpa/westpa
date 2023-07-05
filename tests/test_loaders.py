import os
from unittest import mock
from .test_tools.conftest import *
from westpa.cli.tools.w_multi_west import entry_point, create_idtype_array
from westpa.core.propagators.executable import npy_data_loader, pickle_data_loader
from westpa.core.segment import Segment

import argparse
import numpy as np
import pickle
import h5py


class Test_Loaders:
    '''Class to test w_pdist works to generate a file and that it is the same as the sample pdist.h5 file.'''

    def test_npy_loader(self, ref_idtype):
        '''Test if data loaded with npy_data_loader is consistent.'''

        input_array = h5py.File(self.h5_filepath)['ibstates/0/istate_index'][:]
        test_segment = Segment()

        npy_data_loader('test', self.correct_pkl, test_segment, False) 

        with open(self.correct_pkl, 'rb') as f:
            ref_array = pickle.load(f)

        test_array = test_segment.data['test'][:]

        assert np.array_equal(test_array, ref_array)


    def test_pickle_loader(self, ref_idtype):
        '''Test if data loaded with npy_data_loader is consistent.'''

        input_array = h5py.File(self.h5_filepath)['ibstates/0/istate_index'][:]
        test_segment = Segment()

        pickle_data_loader('test', self.correct_pkl, test_segment, False) 

        with open(self.correct_pkl, 'rb') as f:
            ref_array = pickle.load(f)

        test_array = test_segment.data['test'][:]

        assert np.array_equal(test_array, ref_array)
