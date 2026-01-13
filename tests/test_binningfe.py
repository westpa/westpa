import pytest
from io import StringIO

import numpy as np

from westpa.core.binning.assign import RecursiveBinMapper
from westpa.tools.binning import mapper_from_expr, mapper_from_system, write_bin_info


class TestBinParsing:
    def test_mapper_from_expr(self):
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


class TestBinningWriteout:
    def test_write_bin_info(self):
        out_io = StringIO()
        mapper = mapper_from_expr(r'[[-inf, 0, 15, np.inf]]')
        assignments = mapper.assign([[3], [3], [18]])
        weights = np.asarray([0.25, 0.25, 0.5])

        reference_output = '''3 segments
3 bins total, 1 targets, 2 (100.0%) occupied
Minimum probability by bin:     5.00000000000000000e-01
Maximum probability by bin:     5.00000000000000000e-01
Dynamic range (by bin):         -0 kT
Minimum probability by segment: 2.50000000000000000e-01
Maximum probability by segment: 5.00000000000000000e-01
Dynamic range (by segment):     0.693147 kT
Norm = 1, error in norm = 0 (0 epsilon)

 Index     Count    Total weight               Min seg weight             Max seg weight             Weight ratio    Label
     0         0    0.00000000000000000e+00    0.00000000000000000e+00    0.00000000000000000e+00         0.00000    [(-inf, 0.0)]
     1         2    5.00000000000000000e-01    2.50000000000000000e-01    2.50000000000000000e-01         1.00000    [(0.0, 15.0)]
     2         1    5.00000000000000000e-01    5.00000000000000000e-01    5.00000000000000000e-01         1.00000    [(15.0, inf)]
'''

        write_bin_info(mapper=mapper, assignments=assignments, weights=weights, n_target_states=1, outfile=out_io, detailed=True)
        assert out_io.getvalue() == reference_output
