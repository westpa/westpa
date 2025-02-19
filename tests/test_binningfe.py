import numpy as np

from westpa.tools.binning import mapper_from_expr


class TestBinParsing:
    def testMapperFromExpr(self):
        mapper = mapper_from_expr(r'[[-inf, 0, 15, np.inf]]')

        assert mapper.ndim == 1
        assert mapper.nbins == 3
        assert np.array_equal(mapper._boundaries, [np.array([-float('inf'), 0, 15, float('inf')], dtype=np.float32)])
        assert mapper.labels == ['[(-inf, 0.0)]', '[(0.0, 15.0)]', '[(15.0, inf)]']
