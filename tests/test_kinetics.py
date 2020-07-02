import collections

import numpy as np

from westpa.core.kinetics._kinetics import calc_rates, StreamingStats2D, StreamingStats1D
from westpa.core.kinetics.rate_averaging import tuple2stats


class TestRateAverating:
    def test_tuple2stats2D(self):
        StreamingStatsTuple = collections.namedtuple('StreamingStatsTuple', ['M1', 'M2', 'n'])
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins, nbins))
        mask = np.zeros((nbins, nbins), np.uint8)

        rate_stats = StreamingStats2D((nbins, nbins))
        for d in data:
            rate_stats.update(d, mask)

        c_rate_stats = StreamingStatsTuple(rate_stats.M1, rate_stats.M2, rate_stats.n)

        rate_stats2 = tuple2stats(c_rate_stats)

        assert np.allclose(rate_stats.M1, rate_stats2.M1)
        assert np.allclose(rate_stats.M2, rate_stats2.M2)
        assert np.allclose(rate_stats.n, rate_stats2.n)
        assert np.allclose(rate_stats.mean, rate_stats2.mean)
        assert np.allclose(rate_stats.var, rate_stats2.var)

    def test_tuple2stats1D(self):
        StreamingStatsTuple = collections.namedtuple('StreamingStatsTuple', ['M1', 'M2', 'n'])
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins))
        mask = np.zeros((nbins,), np.uint8)

        rate_stats = StreamingStats1D(nbins)
        for d in data:
            rate_stats.update(d, mask)

        c_rate_stats = StreamingStatsTuple(rate_stats.M1, rate_stats.M2, rate_stats.n)

        rate_stats2 = tuple2stats(c_rate_stats)

        assert np.allclose(rate_stats.M1, rate_stats2.M1)
        assert np.allclose(rate_stats.M2, rate_stats2.M2)
        assert np.allclose(rate_stats.n, rate_stats2.n)
        assert np.allclose(rate_stats.mean, rate_stats2.mean)
        assert np.allclose(rate_stats.var, rate_stats2.var)


class TestKinetics:
    def test_calc_rates(self):
        nbins = 100
        mask = np.zeros((nbins, nbins), np.uint8)

        flux_matrix = np.random.normal(size=(nbins, nbins))
        rate_matrix = np.zeros_like(flux_matrix)

        population_vector = np.random.random(size=(nbins,)) + 0.0001
        population_vector[[0, 2, 5]] = 0.0

        calc_rates(flux_matrix, population_vector, rate_matrix, mask)

        expected_mask = np.zeros((nbins, nbins), np.uint8)
        expected_mask[[0, 2, 5], :] = 1

        nzi = np.where(population_vector != 0.0)[0]
        expected_rate = np.zeros_like(flux_matrix)
        expected_rate[nzi, :] = flux_matrix[nzi, :] / population_vector[nzi][:, np.newaxis]

        assert np.all(mask == expected_mask)
        assert np.allclose(rate_matrix, expected_rate)

        # row_sum = rate_matrix.sum(axis=1)
        # for k in xrange(nbins):
        # print(row_sum[k])
        # assert np.allclose(row_sum[k], 0.0) or np.allclose(row_sum[k], 1.0)


class TestStreamingStats2D:
    def test_nomask(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins, nbins))
        mask = np.zeros((nbins, nbins), np.uint8)

        rate_stats = StreamingStats2D((nbins, nbins))
        for d in data:
            rate_stats.update(d, mask)

        assert np.allclose(rate_stats.mean, data.mean(axis=0))
        assert np.allclose(rate_stats.var, data.var(axis=0,))

    def test_nomask_groups(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins, nbins))
        mask = np.zeros((nbins, nbins), np.uint8)

        rate_stats1 = StreamingStats2D((nbins, nbins))
        for d in data[: (nbins // 2)]:
            rate_stats1.update(d, mask)

        rate_stats2 = StreamingStats2D((nbins, nbins))
        for d in data[(nbins // 2) : nbins]:
            rate_stats2.update(d, mask)

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        assert np.allclose(rate_stats1.mean, data.mean(axis=0))
        assert np.allclose(rate_stats1.var, data.var(axis=0,))

        assert np.allclose(rate_stats3.mean, data.mean(axis=0))
        assert np.allclose(rate_stats3.var, data.var(axis=0,))

    def test_with_mask(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins, nbins))
        mask = np.random.randint(2, size=data.shape).astype(np.uint8)

        rate_stats = StreamingStats2D((nbins, nbins))
        for di, d in enumerate(data):
            rate_stats.update(d, mask[di])

        data_masked = np.ma.array(data, mask=mask)

        assert np.allclose(rate_stats.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert np.allclose(rate_stats.var, data_masked.var(axis=0).filled(fill_value=0.0))

    def test_with_mask_groups(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins, nbins))
        mask = np.random.randint(2, size=data.shape).astype(np.uint8)

        rate_stats1 = StreamingStats2D((nbins, nbins))
        for di, d in enumerate(data[: (nbins // 2)]):
            rate_stats1.update(d, mask[di])

        rate_stats2 = StreamingStats2D((nbins, nbins))
        for di, d in enumerate(data[(nbins // 2) :]):
            rate_stats2.update(d, mask[di])

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        data_masked = np.ma.array(data, mask=mask)

        assert np.allclose(rate_stats1.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert np.allclose(rate_stats1.var, data_masked.var(axis=0).filled(fill_value=0.0))

        assert np.allclose(rate_stats3.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert np.allclose(rate_stats3.var, data_masked.var(axis=0).filled(fill_value=0.0))


class TestStreamingStats1D:
    def test_nomask(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins))
        mask = np.zeros((nbins,), np.uint8)

        rate_stats = StreamingStats1D(nbins)
        for d in data:
            rate_stats.update(d, mask)

        assert np.allclose(rate_stats.mean, data.mean(axis=0))
        assert np.allclose(rate_stats.var, data.var(axis=0,))

    def test_nomask_groups(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins))
        mask = np.zeros((nbins,), np.uint8)

        rate_stats1 = StreamingStats1D(nbins)
        for d in data[: (nbins // 2)]:
            rate_stats1.update(d, mask)

        rate_stats2 = StreamingStats1D(nbins)
        for d in data[(nbins // 2) :]:
            rate_stats2.update(d, mask)

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        assert np.allclose(rate_stats1.mean, data.mean(axis=0))
        assert np.allclose(rate_stats1.var, data.var(axis=0,))

        assert np.allclose(rate_stats3.mean, data.mean(axis=0))
        assert np.allclose(rate_stats3.var, data.var(axis=0,))

    def test_with_mask(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins))
        mask = np.random.randint(2, size=data.shape).astype(np.uint8)

        rate_stats = StreamingStats1D(nbins)
        for di, d in enumerate(data):
            rate_stats.update(d, mask[di])

        data_masked = np.ma.array(data, mask=mask)

        assert np.allclose(rate_stats.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert np.allclose(rate_stats.var, data_masked.var(axis=0).filled(fill_value=0.0))

    def test_with_mask_groups(self):
        nbins = 100
        nsets = 10
        data = np.random.normal(size=(nsets, nbins))
        mask = np.random.randint(2, size=data.shape).astype(np.uint8)

        rate_stats1 = StreamingStats1D(nbins)
        for di, d in enumerate(data[: (nbins // 2)]):
            rate_stats1.update(d, mask[di])

        rate_stats2 = StreamingStats1D(nbins)
        for di, d in enumerate(data[(nbins // 2) :]):
            rate_stats2.update(d, mask[di])

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        data_masked = np.ma.array(data, mask=mask)

        assert np.allclose(rate_stats1.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert np.allclose(rate_stats1.var, data_masked.var(axis=0).filled(fill_value=0.0))

        assert np.allclose(rate_stats3.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert np.allclose(rate_stats3.var, data_masked.var(axis=0).filled(fill_value=0.0))
