import numpy as np

from westpa.core.h5io import WESTIterationFile
from westpa.core.segment import Segment


class Test_H5io:
    '''Class to test the h5io module.'''

    def test_write_segment(self, west_iteration_file):
        '''Test that a segment is written successfully with WestIterationFile.'''

        segment = Segment(n_iter=3, seg_id=5, data=self.dummy_data)
        shape = self.dummy_data['iterh5/trajectory'][:].shape

        with WESTIterationFile(self.h5_iter_file_path, mode='a') as h5_iter_file:
            h5_iter_file.write_segment(segment)
            assert np.allclose(h5_iter_file.root['coordinates'][:], self.dummy_data['iterh5/trajectory'][:])

            self.dummy_data['iterh5/trajectory'][:] *= 2

            # Rewriting with new coordinates and see if it is written in correctly
            h5_iter_file.write_segment(segment)

            assert np.all(h5_iter_file.root['coordinates'][:].shape == shape), 'Extra frames added to per-iter HDF5 File'
            assert np.allclose(
                h5_iter_file.root['coordinates'][:], self.dummy_data['iterh5/trajectory'][:]
            ), 'Unable to overwrite trajectory coordinates'

            # Writing a new segment
            segment2 = Segment(n_iter=3, seg_id=6, data=self.dummy_data)
            h5_iter_file.write_segment(segment2)
            shape = (shape[0] * 2,) + shape[1:]

            assert np.all(
                h5_iter_file.root['coordinates'][:].shape == shape
            ), 'Extra frames not added to per-iter HDF5 File correctly'
            assert np.allclose(
                h5_iter_file.root['coordinates'][2:], self.dummy_data['iterh5/trajectory'][:]
            ), 'Unable to write extra trajectory coordinates'
