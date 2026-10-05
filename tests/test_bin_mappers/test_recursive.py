import numpy as np

from westpa.core.bin_mappers import RecursiveBinMapper
from westpa.core.binning import FuncBinMapper, RectilinearBinMapper
from westpa.core.segment import Segment


class TestRecursiveBinMapper:
    """
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
    """

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
        """
        0                            1                      2
        +----------------------------+----------------------+
        |              0             |          1           |
        +---------------------------------------------------+
        """

        mapper = FuncBinMapper(self.fn1, 2)
        rmapper = RecursiveBinMapper(mapper)
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        segments = [Segment(pcoord=[[0.0], coord]) for coord in coords]
        output = rmapper(segments)
        assert list(output) == [0, 0, 0, 0, 0, 1]

    def testSingleRecursion(self):
        """
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
        """

        outer_mapper = FuncBinMapper(self.fn1, 2)
        inner_mapper = FuncBinMapper(self.fn2, 2)

        rmapper = RecursiveBinMapper(outer_mapper, nested_mappers={0: inner_mapper})

        assert rmapper.nbins == 3
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        segments = [Segment(pcoord=[[0.0], coord]) for coord in coords]
        output = rmapper(segments)
        assert list(output) == [1, 1, 1, 1, 2, 0]

    def testDeepRecursion(self):
        """
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
        """

        outer_mapper = FuncBinMapper(self.fn1, 2)
        middle_mapper = FuncBinMapper(self.fn2, 2)
        inner_mapper = FuncBinMapper(self.fn3, 2)

        rmapper = RecursiveBinMapper(
            outer_mapper,
            nested_mappers={0: RecursiveBinMapper(middle_mapper, nested_mappers={0: inner_mapper})},
        )

        assert rmapper.nbins == 4
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        segments = [Segment(pcoord=[[0.0], coord]) for coord in coords]
        output = rmapper(segments)
        assert list(output) == [2, 2, 3, 3, 1, 0]

    def testSideBySideRecursion(self):
        """
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
        """
        outer_mapper = FuncBinMapper(self.fn1, 2)
        middle_mapper1 = FuncBinMapper(self.fn2, 2)
        middle_mapper2 = FuncBinMapper(self.fn4, 2)

        rmapper = RecursiveBinMapper(
            outer_mapper,
            nested_mappers={0: middle_mapper1, 1: middle_mapper2},
        )

        assert rmapper.nbins == 4
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1], [1.6]])
        segments = [Segment(pcoord=[[0.0], coord]) for coord in coords]
        output = rmapper(segments)
        assert list(output) == [0, 0, 0, 0, 1, 2, 3]

    def testMegaComplexRecursion(self):
        """
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
        """
        outer_mapper = FuncBinMapper(self.fn1, 2)
        middle_mapper1 = FuncBinMapper(self.fn2, 2)
        middle_mapper2 = FuncBinMapper(self.fn4, 2)
        inner_mapper = FuncBinMapper(self.fn3, 2)

        rmapper = RecursiveBinMapper(
            outer_mapper,
            nested_mappers={
                0: RecursiveBinMapper(middle_mapper1, nested_mappers={0: inner_mapper}),
                1: middle_mapper2,
            },
        )

        assert rmapper.nbins == 5
        coords = np.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1], [1.6]])
        segments = [Segment(pcoord=[[0.0], coord]) for coord in coords]
        output = rmapper(segments)
        assert list(output) == [1, 1, 2, 2, 0, 3, 4]

    def test2dRectilinearRecursion(self):
        """
         0                            1                      2
         +----------------------------+----------------------+
         |                            |         1.5          |
         |                            | +--------+---------+ |
         |                            | |        |         | |
         |             0              | |   4    |   5     | |
         |                            | |        |         | |
         |                            | |        |         | |
         |                            | +--------+---------+ |
        1+---------------------------------------------------+
         |            0.5             |                      |
         | +-----------+------------+ |                      |
         | |           |            | |                      |
         | |    2      |     3      | |           1          |
         | |           |            | |                      |
         | |           |            | |                      |
         | +-----------+------------+ |                      |
        2+---------------------------------------------------+

        """
        outer_mapper = RectilinearBinMapper([[0, 1, 2], [0, 1, 2]])

        upper_right_mapper = RectilinearBinMapper([[1, 1.5, 2], [0, 1]])
        lower_left_mapper = RectilinearBinMapper([[0, 0.5, 1], [1, 2]])

        rmapper = RecursiveBinMapper(
            outer_mapper,
            nested_mappers={1: lower_left_mapper, 2: upper_right_mapper},
        )

        pairs = [(0.5, 0.5), (1.25, 0.5), (1.75, 0.5), (0.25, 1.5), (0.75, 1.5), (1.5, 1.5)]
        segments = [Segment(pcoord=[(0.0, 0.0), pair]) for pair in pairs]

        assert rmapper.nbins == 6
        assignments = rmapper(segments)
        labels = list(rmapper.labels)
        expected = [0, 4, 5, 2, 3, 1]
        expected_labels = [
            '[(0.0, 1.0), (0.0, 1.0)]',
            '[(1.0, 2.0), (1.0, 2.0)]',
            '[(0.0, 0.5), (1.0, 2.0)]',
            '[(0.5, 1.0), (1.0, 2.0)]',
            '[(1.0, 1.5), (0.0, 1.0)]',
            '[(1.5, 2.0), (0.0, 1.0)]',
        ]

        assert (assignments == expected).all()
        assert labels == expected_labels
