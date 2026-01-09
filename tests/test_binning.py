import os
import pytest

import h5py
import numpy as np
from scipy.spatial.distance import cdist

import westpa
from westpa.core.binning.assign import (
    RectilinearBinMapper,
    PiecewiseBinMapper,
    FuncBinMapper,
    VectorizingFuncBinMapper,
    VoronoiBinMapper,
    RecursiveBinMapper,
)
from westpa.core.binning.assign import coord_dtype
from westpa.core.binning.mab import MABBinMapper, map_mab, log_bin_boundaries


REFERENCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'refs')


class TestRectilinearBinMapper:
    def test1dAssign(self):
        bounds = [0.0, 1.0, 2.0, 3.0]
        coords = np.array([0, 0.5, 1.5, 1.6, 2.0, 2.0, 2.9])[:, None]

        assigner = RectilinearBinMapper([bounds])
        assert (assigner.assign(coords) == [0, 0, 1, 1, 2, 2, 2]).all()
        assert list(assigner.labels) == ['[(0.0, 1.0)]', '[(1.0, 2.0)]', '[(2.0, 3.0)]']

    def test2dAssign(self):
        """bin structure: [(a,b), (c,d)] => x in [a,b), y in [c, d)"""

        boundaries = [(-1, -0.5, 0, 0.5, 1), (-1, -0.5, 0, 0.5, 1)]
        coords = np.array([(-0.75, -0.75), (-0.25, -0.25), (0, 0), (0.25, 0.25), (0.75, 0.75), (-0.25, 0.75), (0.25, -0.75)])
        assigner = RectilinearBinMapper(boundaries)

        expected_labels = [
            '[(-1.0, -0.5), (-1.0, -0.5)]',
            '[(-1.0, -0.5), (-0.5, 0.0)]',
            '[(-1.0, -0.5), (0.0, 0.5)]',
            '[(-1.0, -0.5), (0.5, 1.0)]',
            '[(-0.5, 0.0), (-1.0, -0.5)]',
            '[(-0.5, 0.0), (-0.5, 0.0)]',
            '[(-0.5, 0.0), (0.0, 0.5)]',
            '[(-0.5, 0.0), (0.5, 1.0)]',
            '[(0.0, 0.5), (-1.0, -0.5)]',
            '[(0.0, 0.5), (-0.5, 0.0)]',
            '[(0.0, 0.5), (0.0, 0.5)]',
            '[(0.0, 0.5), (0.5, 1.0)]',
            '[(0.5, 1.0), (-1.0, -0.5)]',
            '[(0.5, 1.0), (-0.5, 0.0)]',
            '[(0.5, 1.0), (0.0, 0.5)]',
            '[(0.5, 1.0), (0.5, 1.0)]',
        ]

        assert (assigner.assign(coords) == [0, 5, 10, 10, 15, 7, 8]).all()
        assert list(assigner.labels) == expected_labels


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

    def test2dRectilinearRecursion(self):
        '''
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

        print('PAIRS', pairs)
        print('LABELS', labels)
        print('EXPECTED', expected)
        print('OUTPUT  ', assignments)
        assert (assignments == expected).all()
        assert labels == expected_labels


# Following section is for MAB Testing
@pytest.fixture(scope='class')
def gen_input_mab_data_fixture(request):
    '''A fixture to assign `gen_input_mab_data`'s return as a class attribute for tests.'''
    request.cls.input_mab_data = gen_input_mab_data()


def gen_input_mab_data():
    '''Function to generate test data for MABBinMapper Testing.'''
    # Create synthetic test data: a 2D grid of points with Gaussian weights distribution
    N_point_2d = 15
    n_dim_2d = 2
    N_total_2d = N_point_2d**n_dim_2d
    coords_2d = np.zeros((N_total_2d, 2))
    for i in range(N_total_2d):
        coords_2d[i, 0] = i % N_point_2d
        coords_2d[i, 1] = i // N_point_2d
    coords_2d /= N_point_2d

    # generate weights as a 2D gaussian where the max is at the center of the 2D space
    weights_2d = np.zeros(N_total_2d)
    for i in range(N_total_2d):
        weights_2d[i] = np.exp(-((coords_2d[i, 0] - 1 / 2) ** 2 + (coords_2d[i, 1] - 1 / 2) ** 2) / (1 / 2) ** 2)
    weights_2d /= np.sum(weights_2d)

    # Tile coords
    allcoords_2d_grid = np.ones((N_total_2d, 4))
    allcoords_2d_grid[:, :2] = coords_2d
    allcoords_2d_grid[:, 2] = weights_2d
    allcoords_2d_grid = np.tile(allcoords_2d_grid, (2, 1))
    allcoords_2d_grid[0:N_total_2d, 3] = 0

    # Create deterministic 3D array of grid points, again with 3D Gaussian distribution of weights
    N_point_3d = 15
    n_dim_3d = 3
    N_total_3d = N_point_3d**n_dim_3d
    coords_3d = np.zeros((N_total_3d, 3))
    for i in range(N_total_3d):
        coords_3d[i, 0] = i % N_point_3d
        coords_3d[i, 1] = (i // N_point_3d) % N_point_3d
        coords_3d[i, 2] = i // (N_point_3d**2)
    coords_3d /= N_point_3d

    # Generate weights as a n_dim gaussian where the max is at the center of the 3D space
    weights_3d = np.zeros(N_total_3d)
    for i in range(N_total_3d):
        weights_3d[i] = np.exp(
            -((coords_3d[i, 0] - 1 / 2) ** 2 + (coords_3d[i, 1] - 1 / 2) ** 2 + (coords_3d[i, 2] - 1 / 2) ** 2) / (1 / 2) ** 2
        )
    weights_3d /= np.sum(weights_3d)

    # Tile coords
    allcoords_3d_grid = np.ones((N_total_3d, n_dim_3d + 2))
    allcoords_3d_grid[:, :n_dim_3d] = coords_3d
    allcoords_3d_grid[:, n_dim_3d] = weights_3d
    allcoords_3d_grid = np.tile(allcoords_3d_grid, (2, 1))
    allcoords_3d_grid[0:N_total_3d, n_dim_3d + 1] = 0

    # Lastly, test bin assignment with 2D coordinates and weights on pseudorandom deterministic Gaussian points
    N_total_random = 300
    rng = np.random.default_rng(seed=0)
    coords_random = rng.normal(loc=[0.5, 0.5], scale=0.25, size=(N_total_random, n_dim_2d))

    # Generate weights as a n_dim sine curve with given wavelength plus some deterministic noise
    weights_random = np.zeros(N_total_random)
    wavelength = 0.25
    noise_level = 0.1
    for i in range(n_dim_2d):
        weights_random += np.sin(2 * np.pi * coords_random[:, i] / wavelength) + noise_level * np.cos(
            4 * 2 * np.pi * coords_random[:, i] / wavelength
        )
    weights_random = np.abs(weights_random) / np.sum(np.abs(weights_random))

    # Tile coords
    allcoords_2d_gauss = np.ones((N_total_random, 4))
    allcoords_2d_gauss[:, :2] = coords_random
    allcoords_2d_gauss[:, 2] = weights_random
    allcoords_2d_gauss = np.tile(allcoords_2d_gauss, (2, 1))
    allcoords_2d_gauss[0:N_total_random, 3] = 0

    return {
        'allcoords_2d_grid': allcoords_2d_grid,
        'allcoords_3d_grid': allcoords_3d_grid,
        'allcoords_2d_gauss': allcoords_2d_gauss,
    }


@pytest.fixture(scope='class')
def ref_mab_results(request):
    '''Function for reading the reference datasets from `refs/mab_assignments_ref.h5`.'''
    # Setting up empty things dictionary with corresponding keys
    request.cls.ref_mab_results = {}
    test_keys = ['2d_grid', '3d_grid', '2d_gauss']

    # Loading in the three sets of reference data from the reference file into a dictionary
    with h5py.File(os.path.join(REFERENCE_PATH, 'mab_assignments_ref.h5'), 'r') as f:
        for key in test_keys:
            request.cls.ref_mab_results[key] = [f[f'{key}/test_result_{i:d}'][...] for i in range(len(f[key]))]


class TestMABBinMapper:
    '''Test class for MABBinMapper'''

    @pytest.mark.parametrize(
        'nbins, direction, skip, bottleneck, ref_value',
        [
            ([5, 1], [1, 86], [0, 0], True, 9),
            ([5, 5], [0, 0], [0, 0], True, 33),
            ([5, 5, 5], [0, 0, 0], [0, 0, 0], True, 137),
            ([5, 5], [0, 0], [1, 0], True, 9),
            ([5, 1], [1, 86], [0, 0], False, 6),
        ],
        ids=['5x1, direction=[1,86]', '5x5', '5x5x5', '5x5, skip=[1,0]', '5x1, direction=[1,86], no bottleneck'],
    )
    def test_determine_total_bins(self, nbins, direction, skip, bottleneck, ref_value):
        '''Runs through different configurations to see if mab reports the correct number of total bins'''
        mab = MABBinMapper([5])  # Creating a dummy MABBinMapper

        assert mab.determine_total_bins(nbins_per_dim=nbins, direction=direction, skip=skip, bottleneck=bottleneck) == ref_value

    @pytest.mark.parametrize(
        "nbins_per_dim, direction, bottleneck, skip, ref_index",
        [
            ([2, 2], [1, 1], False, [0, 0], 0),
            ([2, 2], [0, 0], False, [0, 0], 1),
            ([2, 2], [-1, -1], True, [0, 0], 2),
            ([2, 2], [0, 0], True, [0, 0], 3),
            ([2, 2], [0, 0], True, [0, 1], 4),
            ([2, 2], [0, 0], True, [1, 1], 5),
            ([2, 2], [86, 0], True, [0, 0], 6),
            ([2, 2], [86, 86], False, [0, 0], 7),
        ],
        ids=[
            'direction=[1,1], no bottleneck',
            'direction=[0,0], no bottleneck',
            'direction=[-1,-1]',
            'direction=[0,0]',
            'direction=[0,0], skip=[0,1]',
            'direction=[0,0], skip=[1,1]',
            'direction=[86,0]',
            'direction=[86,86], no bottleneck',
        ],
    )
    def test_2x2_2d_grid_mab_bin_assignments(
        self, gen_input_mab_data_fixture, ref_mab_results, nbins_per_dim, direction, bottleneck, skip, ref_index
    ):
        '''Test MABBinMapper with 2x2 linear section on 2D space'''
        allcoords = self.input_mab_data['allcoords_2d_grid']
        N_total = allcoords.shape[0] // 2
        mask = np.full((N_total * 2), True)
        output = np.zeros((N_total * 2), dtype=int)
        output = map_mab(
            coords=allcoords,
            mask=mask,
            output=output,
            nbins_per_dim=nbins_per_dim,
            direction=direction,
            bottleneck=bottleneck,
            skip=skip,
        )
        assert np.all(output[:N_total] == output[N_total:]), "Expected first half of bin assignments to equal second half"
        assert np.all(
            output == self.ref_mab_results['2d_grid'][ref_index]
        ), f"Unexpected 2D grid MAB bin assignments with direction={direction}, bottleneck={bottleneck}, and skip={skip}"

    @pytest.mark.parametrize(
        "nbins_per_dim, direction, bottleneck, skip, ref_index",
        [
            ([2, 2, 2], [0, 0, 0], False, [0, 0, 0], 0),
            ([2, 2, 2], [0, 0, 0], True, [0, 0, 0], 1),
        ],
        ids=[
            'no bottleneck',
            'with bottleneck',
        ],
    )
    def test_2x2x2_3d_grid_mab_bin_assignments(
        self, gen_input_mab_data_fixture, ref_mab_results, nbins_per_dim, direction, bottleneck, skip, ref_index
    ):
        '''Test MABBinMapper with 2x2x2 linear section on 3D space'''
        allcoords = self.input_mab_data['allcoords_3d_grid']
        N_total = allcoords.shape[0] // 2
        mask = np.full((N_total * 2), True)
        output = list(np.zeros((N_total * 2), dtype=int))
        output = map_mab(
            coords=allcoords,
            mask=mask,
            output=output,
            nbins_per_dim=nbins_per_dim,
            direction=direction,
            bottleneck=bottleneck,
            skip=skip,
        )
        assert output[:N_total] == output[N_total:], "Expected first half of bin assignments to equal second half"
        assert output == list(
            self.ref_mab_results['3d_grid'][ref_index]
        ), f"Unexpected 3D grid MAB bin assignments with direction={direction}, bottleneck={bottleneck}, and skip={skip}"

    @pytest.mark.parametrize(
        "nbins_per_dim, direction, bottleneck, skip, ref_index",
        [
            ([2, 2], [0, 0], True, [0, 0], 0),
            ([2, 2], [0, 0], True, [0, 1], 1),
            ([2, 2], [86, -1], True, [0, 0], 2),
        ],
        ids=[
            'direction=[0,0], no skip',
            'direction=[0,0], skip=[0,1]',
            'direction=[86,-1], no skip',
        ],
    )
    def test_2x2_2d_gaussian_mab_bin_assignments(
        self, gen_input_mab_data_fixture, ref_mab_results, nbins_per_dim, direction, bottleneck, skip, ref_index
    ):
        '''Test MABBinMapper with 2x2 linear section on a 2D Gaussian space'''
        allcoords = self.input_mab_data['allcoords_2d_gauss']
        N_total = allcoords.shape[0] // 2
        mask = np.full((N_total * 2), True)
        output = np.zeros((N_total * 2), dtype=int)
        output = map_mab(
            coords=allcoords,
            mask=mask,
            output=output,
            nbins_per_dim=nbins_per_dim,
            direction=direction,
            bottleneck=bottleneck,
            skip=skip,
        )
        assert np.all(output[:N_total] == output[N_total:]), "Expected first half of bin assignments to equal second half"
        assert np.all(
            output == self.ref_mab_results['2d_gauss'][ref_index]
        ), f"Unexpected 2D Gaussian MAB bin assignments with direction={direction}, bottleneck={bottleneck}, and skip={skip}"

    @pytest.mark.parametrize(
        "skip, bottleneck, direction, minlist, maxlist, nbins_per_dim, n_bottleneck_filled, bottlenecks_forward, bottlenecks_reverse",
        [
            ([0, 0], True, [0, -1], [0.0, 0.0], [1.0, 1.0], [2, 2], 0, [None, None], [None, None]),
            ([0, 0], True, [0, -1], [0.0, 0.0], [1.0, 1.0], [2, 2], 2, [1.0, 1.0], [0.0, 0.0]),
            ([0, 0], False, [0, -1], [0.0, 0.0], [1.0, 1.0], [2, 2], 2, [], []),
        ],
        ids=[
            'None as bottlenecks',
            'With bottlenecks',
            'No bottlenecks',
        ],
    )
    def test_log_bin_boundaries(
        self,
        ref_mab,
        monkeypatch,
        skip,
        bottleneck,
        direction,
        minlist,
        maxlist,
        nbins_per_dim,
        n_bottleneck_filled,
        bottlenecks_forward,
        bottlenecks_reverse,
    ):
        '''Test MAB logging with various situations'''

        temp_file_path = f'{self.tmpdir}/log_output.txt'

        with monkeypatch.context() as m:
            m.setattr(westpa, 'rc', westpa.core._rc.WESTRC())
            westpa.rc.read_config(filename='west.cfg')
            m.setattr(westpa.rc.sim_manager, 'n_iter', 24)

            log_bin_boundaries(
                skip=skip,
                bottleneck=bottleneck,
                direction=direction,
                bin_log_path=temp_file_path,
                minlist=minlist,
                maxlist=maxlist,
                nbins_per_dim=nbins_per_dim,
                n_bottleneck_filled=n_bottleneck_filled,
                bottlenecks_forward=bottlenecks_forward,
                bottlenecks_reverse=bottlenecks_reverse,
            )

        # Correct outputs for comparison.
        template_output = '''Iteration: 24
MAB linear bin boundaries: [0.  0.5 1. ]\t[0.  0.5 1. ]\t
Lagging pcoord in each dimension: [0.0, 0.0]
Leading pcoord in each dimension: [1.0, 1.0]
'''

        if bottleneck:
            template_output += f'''Number of bottleneck bins filled: {n_bottleneck_filled} / 3
Dimension 0 forward bottleneck walker at: [{bottlenecks_forward[0]}]
Dimension 0 backward bottleneck walker at: [{bottlenecks_reverse[0]}]
Dimension 1 backward bottleneck walker at: [{bottlenecks_reverse[1]}]

'''
        else:
            template_output += '\n'

        # Do the actual comparison...
        with open(temp_file_path, 'r') as f:
            assert template_output == f.read()


def output_mab_reference():
    '''
    Function to generate the reference test files for the MAB tests.

    To run this, run `python -c 'from test_binning import output_mab_reference; output_mab_reference()'`
    in the command line (assuming you're in the `tests/` folder).

    It will overwrite the file 'tests/refs/mab_assignments_ref.h5' and create 11 png files in `tests/` visualizing the binning.
    '''
    import matplotlib.pyplot as plt

    try:
        os.mkdir('mab_tests')
    except FileExistsError:
        pass

    input_data = gen_input_mab_data()
    with h5py.File(f'{REFERENCE_PATH}/mab_assignments_ref.h5', 'w') as f:
        # 2D Grid
        for i, (nbins_per_dim, direction, bottleneck, skip) in enumerate(
            [
                ([2, 2], [1, 1], False, [0, 0]),
                ([2, 2], [0, 0], False, [0, 0]),
                ([2, 2], [-1, -1], True, [0, 0]),
                ([2, 2], [0, 0], True, [0, 0]),
                ([2, 2], [0, 0], True, [0, 1]),
                ([2, 2], [0, 0], True, [1, 1]),
                ([2, 2], [86, 0], True, [0, 0]),
                ([2, 2], [86, 86], False, [0, 0]),
            ]
        ):
            allcoords = input_data['allcoords_2d_grid']
            N_total = allcoords.shape[0] // 2
            mask = np.full((N_total * 2), True)
            output = np.zeros((N_total * 2), dtype=int)
            output = map_mab(
                coords=allcoords,
                mask=mask,
                output=output,
                nbins_per_dim=nbins_per_dim,
                direction=direction,
                bottleneck=bottleneck,
                skip=skip,
            )

            f.create_dataset(f'2d_grid/test_result_{i:d}', data=output)

            # Create a cmap with the same number of colors as the number of bins
            cmap = plt.cm.get_cmap('tab20', int(np.max(output) + 1))

            # Plot the synthetic data in 2D using a scatter plot
            # Include a cbar to shown the bin assignments
            plt.scatter(
                allcoords[:N_total, 0],
                allcoords[:N_total, 1],
                s=allcoords[:N_total, 2] * 10000,
                c=output[:N_total],
                cmap=cmap,
                vmin=-0.5,
                vmax=int(np.max(output)) + 0.5,
            )
            plt.colorbar(label='Bin ID', ticks=range(int(np.max(output) + 1)))
            plt.title(f'nbins_per_dim={nbins_per_dim}, direction={direction},\n bottleneck={bottleneck}, skip={skip}')
            plt.savefig(f'mab_tests/2d_grid_ref_result_{i}.png')
            plt.clf()

        # 3D grid
        for i, (nbins_per_dim, direction, bottleneck, skip) in enumerate(
            [
                ([2, 2, 2], [0, 0, 0], False, [0, 0, 0]),
                ([2, 2, 2], [0, 0, 0], True, [0, 0, 0]),
            ]
        ):
            allcoords = input_data['allcoords_3d_grid']
            N_total = allcoords.shape[0] // 2
            mask = np.full((N_total * 2), True)
            output = np.zeros((N_total * 2), dtype=int)
            output = map_mab(
                coords=allcoords,
                mask=mask,
                output=output,
                nbins_per_dim=nbins_per_dim,
                direction=direction,
                bottleneck=bottleneck,
                skip=skip,
            )

            f.create_dataset(f'3d_grid/test_result_{i:d}', data=output)

        # 2D Gaussian
        for i, (nbins_per_dim, direction, bottleneck, skip) in enumerate(
            [
                ([2, 2], [0, 0], True, [0, 0]),
                ([2, 2], [0, 0], True, [0, 1]),
                ([2, 2], [86, -1], True, [0, 0]),
            ]
        ):
            allcoords = input_data['allcoords_2d_gauss']
            N_total = allcoords.shape[0] // 2
            mask = np.full((N_total * 2), True)
            output = np.zeros((N_total * 2), dtype=int)
            output = map_mab(
                coords=allcoords,
                mask=mask,
                output=output,
                nbins_per_dim=nbins_per_dim,
                direction=direction,
                bottleneck=bottleneck,
                skip=skip,
            )

            f.create_dataset(f'2d_gauss/test_result_{i:d}', data=output)

            # Create a cmap with the same number of colors as the number of bins
            cmap = plt.cm.get_cmap('tab20', int(np.max(output) + 1))

            # Plot the synthetic data in 2D using a scatter plot
            # Include a cbar to shown the bin assignments
            plt.scatter(
                allcoords[:N_total, 0],
                allcoords[:N_total, 1],
                s=allcoords[:N_total, 2] * 10000,
                c=output[:N_total],
                cmap=cmap,
                vmin=-0.5,
                vmax=int(np.max(output)) + 0.5,
            )
            plt.colorbar(label='Bin ID', ticks=range(int(np.max(output) + 1)))
            plt.title(f'nbins_per_dim={nbins_per_dim}, direction={direction},\n bottleneck={bottleneck}, skip={skip}')
            plt.savefig(f'mab_tests/2d_gauss_ref_result_{i}.png')
            plt.clf()
    print("Reference data generated and saved to file.")
