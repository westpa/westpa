import numpy as np

from westpa.core.bin_mappers import BinMapperBase
from westpa.core.segment import Segment


class TestBinMapperBase:

    def test_two_bin_implementation(self):
        class TwoBinMapper(BinMapperBase):
            @property
            def labels(self):
                return ['-', '+']

            def map(self, coords, weights, output):
                output[:] = coords[:, 0] > 0
                return output

        bin_mapper = TwoBinMapper()

        assert bin_mapper.labels == ['-', '+']
        assert bin_mapper.nbins == 2

        pcoord = np.array([[-1.0], [0.0], [1.0]])
        segments = [Segment(pcoord=pcoord), Segment(pcoord=-pcoord)]
        assert bin_mapper(segments).tolist() == [1, 0]
        assert bin_mapper(segments, coord_index=0).tolist() == [0, 1]
        assert bin_mapper(segments, coord_index=1).tolist() == [0, 0]
