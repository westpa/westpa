# Copyright (C) 2013 Joshua L. Adelman
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.


from westpa.kinetics._kinetics import calc_rates, StreamingStats2D, StreamingStats1D
from westpa.kinetics.rate_averaging import tuple2stats
from collections import namedtuple

import numpy
import nose
import nose.tools


class TestRateAverating:
    def test_tuple2stats2D(self):
        StreamingStatsTuple = namedtuple('StreamingStatsTuple', ['M1', 'M2', 'n'])
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        mask = numpy.zeros((nbins, nbins), numpy.uint8)

        rate_stats = StreamingStats2D((nbins, nbins))
        for d in data:
            rate_stats.update(d, mask)

        c_rate_stats = StreamingStatsTuple(rate_stats.M1, rate_stats.M2, rate_stats.n)

        rate_stats2 = tuple2stats(c_rate_stats)

        assert numpy.allclose(rate_stats.M1, rate_stats2.M1)
        assert numpy.allclose(rate_stats.M2, rate_stats2.M2)
        assert numpy.allclose(rate_stats.n, rate_stats2.n)
        assert numpy.allclose(rate_stats.mean, rate_stats2.mean)
        assert numpy.allclose(rate_stats.var, rate_stats2.var)

    def test_tuple2stats1D(self):
        StreamingStatsTuple = namedtuple('StreamingStatsTuple', ['M1', 'M2', 'n'])
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        mask = numpy.zeros((nbins,), numpy.uint8)

        rate_stats = StreamingStats1D(nbins)
        for d in data:
            rate_stats.update(d, mask)

        c_rate_stats = StreamingStatsTuple(rate_stats.M1, rate_stats.M2, rate_stats.n)

        rate_stats2 = tuple2stats(c_rate_stats)

        assert numpy.allclose(rate_stats.M1, rate_stats2.M1)
        assert numpy.allclose(rate_stats.M2, rate_stats2.M2)
        assert numpy.allclose(rate_stats.n, rate_stats2.n)
        assert numpy.allclose(rate_stats.mean, rate_stats2.mean)
        assert numpy.allclose(rate_stats.var, rate_stats2.var)

class TestKinetics:
    def test_calc_rates(self):
        nbins = 100
        mask = numpy.zeros((nbins,nbins), numpy.uint8)

        flux_matrix = numpy.random.normal(size=(nbins,nbins))
        rate_matrix = numpy.zeros_like(flux_matrix)

        population_vector = numpy.random.random(size=(nbins,)) + 0.0001
        population_vector[[0,2,5]] = 0.0

        calc_rates(flux_matrix, population_vector, rate_matrix, mask)

        expected_mask = numpy.zeros((nbins,nbins), numpy.uint8)
        expected_mask[[0,2,5],:] = 1

        nzi = numpy.where(population_vector != 0.0)[0]
        expected_rate = numpy.zeros_like(flux_matrix)
        expected_rate[nzi,:] = flux_matrix[nzi,:] / population_vector[nzi][:,numpy.newaxis]

        assert numpy.all(mask == expected_mask)
        assert numpy.allclose(rate_matrix, expected_rate)

        #row_sum = rate_matrix.sum(axis=1)
        #for k in xrange(nbins):
            #print(row_sum[k])
            #assert numpy.allclose(row_sum[k], 0.0) or numpy.allclose(row_sum[k], 1.0)

class TestStreamingStats2D:
    def test_nomask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        mask = numpy.zeros((nbins, nbins), numpy.uint8)

        rate_stats = StreamingStats2D((nbins, nbins))
        for d in data:
            rate_stats.update(d, mask)

        assert numpy.allclose(rate_stats.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats.var, data.var(axis=0,))

    def test_nomask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        mask = numpy.zeros((nbins, nbins), numpy.uint8)

        rate_stats1 = StreamingStats2D((nbins, nbins))
        for d in data[:(nbins/2)]:
            rate_stats1.update(d, mask)

        rate_stats2 = StreamingStats2D((nbins, nbins))
        for d in data[(nbins/2):nbins]:
            rate_stats2.update(d, mask)

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        assert numpy.allclose(rate_stats1.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats1.var, data.var(axis=0,))

        assert numpy.allclose(rate_stats3.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats3.var, data.var(axis=0,))

    def test_with_mask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats = StreamingStats2D((nbins, nbins))
        for di, d in enumerate(data):
            rate_stats.update(d, mask[di])

        data_masked = numpy.ma.array(data, mask=mask)

        assert numpy.allclose(rate_stats.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert numpy.allclose(rate_stats.var, data_masked.var(axis=0).filled(fill_value=0.0))

    def test_with_mask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats1 = StreamingStats2D((nbins, nbins))
        for di, d in enumerate(data[:(nbins/2)]):
            rate_stats1.update(d, mask[di])

        rate_stats2 = StreamingStats2D((nbins, nbins))
        for di, d in enumerate(data[(nbins/2):]):
            rate_stats2.update(d, mask[di])

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        data_masked = numpy.ma.array(data, mask=mask)

        assert numpy.allclose(rate_stats1.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert numpy.allclose(rate_stats1.var, data_masked.var(axis=0).filled(fill_value=0.0))

        assert numpy.allclose(rate_stats3.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert numpy.allclose(rate_stats3.var, data_masked.var(axis=0).filled(fill_value=0.0))


class TestStreamingStats1D:
    def test_nomask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        mask = numpy.zeros((nbins,), numpy.uint8)

        rate_stats = StreamingStats1D(nbins)
        for d in data:
            rate_stats.update(d, mask)

        assert numpy.allclose(rate_stats.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats.var, data.var(axis=0,))

    def test_nomask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        mask = numpy.zeros((nbins,), numpy.uint8)

        rate_stats1 = StreamingStats1D(nbins)
        for d in data[:(nbins/2)]:
            rate_stats1.update(d, mask)

        rate_stats2 = StreamingStats1D(nbins)
        for d in data[(nbins/2):]:
            rate_stats2.update(d, mask)

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        assert numpy.allclose(rate_stats1.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats1.var, data.var(axis=0,))

        assert numpy.allclose(rate_stats3.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats3.var, data.var(axis=0,))

    def test_with_mask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats = StreamingStats1D(nbins)
        for di, d in enumerate(data):
            rate_stats.update(d, mask[di])

        data_masked = numpy.ma.array(data, mask=mask)

        assert numpy.allclose(rate_stats.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert numpy.allclose(rate_stats.var, data_masked.var(axis=0).filled(fill_value=0.0))

    def test_with_mask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats1 = StreamingStats1D(nbins)
        for di, d in enumerate(data[:(nbins/2)]):
            rate_stats1.update(d, mask[di])

        rate_stats2 = StreamingStats1D(nbins)
        for di, d in enumerate(data[(nbins/2):]):
            rate_stats2.update(d, mask[di])

        rate_stats3 = rate_stats1 + rate_stats2
        rate_stats1 += rate_stats2

        data_masked = numpy.ma.array(data, mask=mask)

        assert numpy.allclose(rate_stats1.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert numpy.allclose(rate_stats1.var, data_masked.var(axis=0).filled(fill_value=0.0))

        assert numpy.allclose(rate_stats3.mean, data_masked.mean(axis=0).filled(fill_value=0.0))
        assert numpy.allclose(rate_stats3.var, data_masked.var(axis=0).filled(fill_value=0.0))
