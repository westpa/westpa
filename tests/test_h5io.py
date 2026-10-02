import numpy as np

from westpa.core.h5io import WESTIterationFile
from westpa.core.segment import Segment


class Test_H5io:
    """Class to test the h5io module."""

    def test_write_segment(self, west_iteration_file):
        """Segment is written successfully with WestIterationFile when there are duplicate pointers (and not)."""

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
            old_shape0 = shape[0]
            shape = (shape[0] * 2,) + shape[1:]

            assert np.all(
                h5_iter_file.root['coordinates'][:].shape == shape
            ), 'Extra frames not added to per-iter HDF5 File correctly'
            assert np.allclose(
                h5_iter_file.root['coordinates'][-old_shape0:], self.dummy_data['iterh5/trajectory'][:]
            ), 'Unable to write extra trajectory coordinates'

    def test_write_segment_under(self, west_iteration_file):
        """Segment is written successfully in WestIterationFile when there are not enough rows."""

        oneframe_data = {'iterh5/trajectory': self.dummy_data['iterh5/trajectory'][:1].copy()}
        self.dummy_data['iterh5/trajectory'][:] *= 2

        segment = Segment(n_iter=3, seg_id=5, data=self.dummy_data)
        segment_oneframe = Segment(n_iter=3, seg_id=5, data=oneframe_data)
        shape = self.dummy_data['iterh5/trajectory'][:].shape

        with WESTIterationFile(self.h5_iter_file_path, mode='a') as h5_iter_file:
            h5_iter_file.write_segment(segment_oneframe)
            assert np.allclose(h5_iter_file.root['coordinates'][:], oneframe_data['iterh5/trajectory'][:])
            print(self.dummy_data['iterh5/trajectory'][:])
            print(h5_iter_file.root['coordinates'][:])
            print(h5_iter_file.root['pointer'][:])
            # Rewriting with new coordinates and see if it is written in correctly
            h5_iter_file.write_segment(segment)

            assert np.all(h5_iter_file.root['coordinates'][:].shape == shape), 'Extra frames added to per-iter HDF5 File'
            assert np.allclose(
                h5_iter_file.root['coordinates'][:], self.dummy_data['iterh5/trajectory'][:]
            ), 'Unable to overwrite trajectory coordinates'

    def test_write_segment_over(self, west_iteration_file):
        """Segment is written successfully in WestIterationFile when there are extra rows."""

        oneframe_data = {'iterh5/trajectory': self.dummy_data['iterh5/trajectory'][:1].copy()}
        self.dummy_data['iterh5/trajectory'][:] *= 2

        segment = Segment(n_iter=3, seg_id=5, data=self.dummy_data)
        segment_oneframe = Segment(n_iter=3, seg_id=5, data=oneframe_data)
        shape = self.dummy_data['iterh5/trajectory'][:].shape

        with WESTIterationFile(self.h5_iter_file_path, mode='a') as h5_iter_file:
            h5_iter_file.write_segment(segment)
            assert np.allclose(h5_iter_file.root['coordinates'][:], self.dummy_data['iterh5/trajectory'][:])
            print(self.dummy_data['iterh5/trajectory'][:])
            print(h5_iter_file.root['coordinates'][:])
            # Rewriting with new coordinates and see if it is written in correctly
            h5_iter_file.write_segment(segment_oneframe)

            assert np.all(h5_iter_file.root['coordinates'][:].shape == shape), 'Extra frames added to per-iter HDF5 File'
            assert np.allclose(
                h5_iter_file.root['coordinates'][:1], oneframe_data['iterh5/trajectory'][:]
            ), 'Unable to overwrite trajectory coordinates'
