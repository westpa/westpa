import pytest

import numpy as np

from westpa.tools.binning import mapper_from_expr, mapper_from_system
from westpa.core.binning.assign import RecursiveBinMapper


class TestBinParsing:
    def testMapperFromExpr(self):
        mapper = mapper_from_expr(r'[[-inf, 0, 15, np.inf]]')

        assert mapper.ndim == 1
        assert mapper.nbins == 3
        assert np.array_equal(mapper._boundaries, [np.array([-float('inf'), 0, 15, float('inf')], dtype=np.float32)])
        assert mapper.labels == ['[(-inf, 0.0)]', '[(0.0, 15.0)]', '[(15.0, inf)]']

    def test_mapper_from_expr_exceptions(self):
        with pytest.raises(NameError):
            mapper_from_expr('abc')

        with pytest.raises(TypeError):
            mapper_from_expr(int(3))

    def test_mapper_from_system(self, ref_mab):
        assert isinstance(mapper_from_system(), RecursiveBinMapper)
