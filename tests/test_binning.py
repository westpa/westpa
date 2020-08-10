import pytest

import numpy as np
from scipy.spatial.distance import cdist

from westpa.core.binning.assign import (
    RectilinearBinMapper,
    PiecewiseBinMapper,
    FuncBinMapper,
    VectorizingFuncBinMapper,
    VoronoiBinMapper,
    RecursiveBinMapper,
)
from westpa.core.binning.assign import coord_dtype


class TestRectilinearBinMapper:
    def test1dAssign(self):
        bounds = [0.0, 1.0, 2.0, 3.0]
        coords = np.array([0, 0.5, 1.5, 1.6, 2.0, 2.0, 2.9])[:, None]

        assigner = RectilinearBinMapper([bounds])
        assert (assigner.assign(coords) == [0, 0, 1, 1, 2, 2, 2]).all()

    def test2dAssign(self):
        boundaries = [(-1, -0.5, 0, 0.5, 1), (-1, -0.5, 0, 0.5, 1)]
        coords = np.array([(-0.75, -0.75), (-0.25, -0.25), (0, 0), (0.25, 0.25), (0.75, 0.75), (-0.25, 0.75), (0.25, -0.75)])
        assigner = RectilinearBinMapper(boundaries)

        """bin structure: [(a,b), (c,d)] => x in [a,b), y in [c, d)
        0:[(-1, -0.5), (-1, -0.5)]
        1:[(-1, -0.5), (-0.5, 0)]
        2:[(-1, -0.5), (0, 0.5)]
        3:[(-1, -0.5), (0.5, 1)]
        4:[(-0.5, 0), (-1, -0.5)]
        5:[(-0.5, 0), (-0.5, 0)]
        6:[(-0.5, 0), (0, 0.5)]
        7:[(-0.5, 0), (0.5, 1)]
        8:[(0, 0.5), (-1, -0.5)]
        9:[(0, 0.5), (-0.5, 0)]
        10:[(0, 0.5), (0, 0.5)]
        11:[(0, 0.5), (0.5, 1)]
        12:[(0.5, 1), (-1, -0.5)]
        13:[(0.5, 1), (-0.5, 0)]
        14:[(0.5, 1), (0, 0.5)]
        15:[(0.5, 1), (0.5, 1)]"""

        assert (assigner.assign(coords) == [0, 5, 10, 10, 15, 7, 8]).all()


class TestPiecewiseBinMapper:
    def test_bin_mapping(self):
        coords = np.array([[-0.5], [0.0], [0.5]], dtype=np.float32)
        fr1 = lambda x: x < 0
        fr2 = lambda x: x >= 0
        pm = PiecewiseBinMapper([fr1, fr2])
        assert list(pm.assign(coords)) == [0, 1, 1]


class TestFuncBinMapper:
    @staticmethod
    def fn(coords, mask, output):
        output[mask & (coords[:, 0] <= 0.5)] = 0
        output[mask & (coords[:, 0] > 0.5)] = 1

    def test_fmapper(self):
        mapper = FuncBinMapper(self.fn, 2)
        coords = np.array([0.0, 0.1, 0.5, 0.7])
        coords.shape = (coords.shape[0], 1)
        print(repr(coords))
        output = mapper.assign(coords)
        print(repr(output))
        assert list(output) == [0, 0, 0, 1]


class TestVectorizingFuncBinMapper:
    @staticmethod
    def fn(coord):
        if coord[0] <= 0.5:
            return 0
        else:
            return 1

    def test_vmapper(self):
        mapper = VectorizingFuncBinMapper(self.fn, 2)
        coords = np.array([0.0, 0.1, 0.5, 0.7])
        coords.shape = (coords.shape[0], 1)
        output = mapper.assign(coords)
        assert list(output) == [0, 0, 0, 1]


class TestVoronoiBinMapper:
    @staticmethod
    def distfunc(coordvec, centers):
        if coordvec.ndim < 2:
            new_coordvec = np.empty((1, coordvec.shape[0]), dtype=coord_dtype)
            new_coordvec[0, :] = coordvec[:]
            coordvec = new_coordvec
        distmat = np.require(cdist(coordvec, centers), dtype=coord_dtype)
        return distmat[0, :]

    def test_vmapper(self):
        centers = np.array([[0, 0], [2, 2]], dtype=coord_dtype)
        coords = np.array([[0, 0], [2, 2], [0.9, 0.9], [1.1, 1.1]], dtype=coord_dtype)

        mapper = VoronoiBinMapper(self.distfunc, centers)
        output = mapper.assign(coords)
        assert list(output) == [0, 1, 0, 1]


class TestNestingBinMapper:
    # pass
    '''
         0                            1                      2
         +----------------------------+----------------------+
         |            0.5             |         1.5          |
         | +-----------+------------+ | +--------+---------+ |
         | |    0.25   |            | | |        |         | |
         | | +---+---+ |            | | |        |         | |
         | | |   |   | |            | | |        |         | |
         | | +---+---+ |            | | |        |         | |
         | +-----------+------------+ | +--------+---------+ |
         +---------------------------------------------------+
    '''

    @staticmethod
    def fn1(coords, mask, output):
        test = coords[:, 0] < 1
        output[mask & test] = 0
        output[mask & ~test] = 1

    @staticmethod
    def fn2(coords, mask, output):
        test = coords[:, 0] < 0.5
        output[mask & test] = 0
        output[mask & ~test] = 1

    @staticmethod
    def fn3(coords, mask, output):
        test = coords[:, 0] < 0.25
        output[mask & test] = 0
        output[mask & ~test] = 1

    @staticmethod
    def fn4(coords, mask, output):
        test = coords[:, 0] < 1.5
        output[mask & test] = 0
        output[mask & ~test] = 1

    def testOuterMapper(self):
        '''
             0                            1                      2
             +----------------------------+----------------------+
             |              0             |          1           |
             +---------------------------------------------------+
        '''

        mapper = FuncBinMapper(self.fn1, 2)
        rmapper = RecursiveBinMapper(mapper)
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [0, 0, 0, 0, 0, 1]

    def testOuterMapperWithOffset(self):
        '''
             0                            1                      2
             +----------------------------+----------------------+
             |              1             |          2           |
             +---------------------------------------------------+
        '''

        mapper = FuncBinMapper(self.fn1, 2)
        rmapper = RecursiveBinMapper(mapper, 1)
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [1, 1, 1, 1, 1, 2]

    def testSingleRecursion(self):
        '''
             0                            1                      2
             +----------------------------+----------------------+
             |            0.5             |                      |
             | +-----------+------------+ |                      |
             | |           |            | |                      |
             | |     1     |     2      | |          0           |
             | |           |            | |                      |
             | |           |            | |                      |
             | +-----------+------------+ |                      |
             +---------------------------------------------------+
        '''

        outer_mapper = FuncBinMapper(self.fn1, 2)
        inner_mapper = FuncBinMapper(self.fn2, 2)
        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(inner_mapper, [0.5])
        assert rmapper.nbins == 3
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [1, 1, 1, 1, 2, 0]

    def testDeepRecursion(self):
        '''
         0                            1                      2
         +----------------------------+----------------------+
         |            0.5             |                      |
         | +-----------+------------+ |                      |
         | |    0.25   |            | |                      |
         | | +---+---+ |            | |           0          |
         | | | 2 | 3 | |     1      | |                      |
         | | +---+---+ |            | |                      |
         | +-----------+------------+ |                      |
         +---------------------------------------------------+
        '''

        outer_mapper = FuncBinMapper(self.fn1, 2)
        middle_mapper = FuncBinMapper(self.fn2, 2)
        inner_mapper = FuncBinMapper(self.fn3, 2)

        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(middle_mapper, [0.5])
        rmapper.add_mapper(inner_mapper, [0.25])

        assert rmapper.nbins == 4
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [2, 2, 3, 3, 1, 0]

    def testSideBySideRecursion(self):
        '''
             0                            1                      2
             +----------------------------+----------------------+
             |            0.5             |         1.5          |
             | +-----------+------------+ | +--------+---------+ |
             | |           |            | | |        |         | |
             | |     0     |     1      | | |   2    |    3    | |
             | |           |            | | |        |         | |
             | |           |            | | |        |         | |
             | +-----------+------------+ | +--------+---------+ |
             +---------------------------------------------------+
        '''
        outer_mapper = FuncBinMapper(self.fn1, 2)
        middle_mapper1 = FuncBinMapper(self.fn2, 2)
        middle_mapper2 = FuncBinMapper(self.fn4, 2)

        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(middle_mapper1, [0.5])
        rmapper.add_mapper(middle_mapper2, [1.5])

        assert rmapper.nbins == 4
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1], [1.6]])
        output = rmapper.assign(coords)
        print('OUTPUT', output)
        assert list(output) == [0, 0, 0, 0, 1, 2, 3]

    def testMegaComplexRecursion(self):
        '''
             0                            1                      2
             +----------------------------+----------------------+
             |            0.5             |         1.5          |
             | +-----------+------------+ | +--------+---------+ |
             | |    0.25   |            | | |        |         | |
             | | +---+---+ |            | | |        |         | |
             | | | 1 | 2 | |     0      | | |   3    |    4    | |
             | | +---+---+ |            | | |        |         | |
             | +-----------+------------+ | +--------+---------+ |
             +---------------------------------------------------+
        '''
        outer_mapper = FuncBinMapper(self.fn1, 2)
        middle_mapper1 = FuncBinMapper(self.fn2, 2)
        middle_mapper2 = FuncBinMapper(self.fn4, 2)
        inner_mapper = FuncBinMapper(self.fn3, 2)

        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(middle_mapper1, [0.5])
        rmapper.add_mapper(inner_mapper, [0.25])
        rmapper.add_mapper(middle_mapper2, [1.5])

        assert rmapper.nbins == 5
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1], [1.6]])
        output = rmapper.assign(coords)
        assert list(output) == [1, 1, 2, 2, 0, 3, 4]

    # TODO: Fix this test
    @pytest.mark.xfail(reason="known error in assign")
    def test2dRectilinearRecursion(self):

        '''
             0                            1                      2
             +----------------------------+----------------------+
             |                            |         1.5          |
             |                            | +--------+---------+ |
             |                            | |        |         | |
             |             0              | |   3    |   4     | |
             |                            | |        |         | |
             |                            | |        |         | |
             |                            | +--------+---------+ |
            1+---------------------------------------------------+
             |            0.5             |                      |
             | +-----------+------------+ |                      |
             | |           |            | |                      |
             | |    1      |     2      | |           5          |
             | |           |            | |                      |
             | |           |            | |                      |
             | +-----------+------------+ |                      |
            2+---------------------------------------------------+

        '''

        outer_mapper = RectilinearBinMapper([[0, 1, 2], [0, 1, 2]])

        upper_right_mapper = RectilinearBinMapper([[1, 1.5, 2], [0, 1]])
        lower_left_mapper = RectilinearBinMapper([[0, 0.5, 1], [1, 2]])

        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(upper_right_mapper, [1.5, 0.5])
        rmapper.add_mapper(lower_left_mapper, [0.5, 1.5])

        pairs = [(0.5, 0.5), (1.25, 0.5), (1.75, 0.5), (0.25, 1.5), (0.75, 1.5), (1.5, 1.5)]

        assert rmapper.nbins == 6
        assignments = rmapper.assign(pairs)
        expected = [0, 3, 4, 1, 2, 5]
        print('PAIRS', pairs)
        print('LABELS', list(rmapper.labels))
        print('EXPECTED', expected)
        print('OUTPUT  ', assignments)
        assert (assignments == expected).all()
