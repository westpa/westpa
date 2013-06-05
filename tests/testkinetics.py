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

from __future__ import division, print_function
from westpa.kinetics._kinetics import calc_rates, StreamingStats2D, StreamingStats1D

import numpy
import nose
import nose.tools


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

class TestStreamingStats2D:
    def test_nomask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        n_vals = numpy.ones((nbins, nbins), numpy.uint)
        mask = numpy.zeros((nbins, nbins), numpy.uint8)

        rate_stats = StreamingStats2D(nbins)
        for d in data:
            rate_stats.update(d, d**2, n_vals, mask)

        assert numpy.allclose(rate_stats.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats.var, data.var(axis=0,))

    def test_nomask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        n_vals = numpy.ones((nbins, nbins), numpy.uint)
        mask = numpy.zeros((nbins, nbins), numpy.uint8)

        rate_stats1 = StreamingStats2D(nbins)
        for d in data[:50]:
            rate_stats1.update(d, d**2, n_vals, mask)

        rate_stats2 = StreamingStats2D(nbins)
        for d in data[50:]:
            rate_stats2.update(d, d**2, n_vals, mask)

        rate_stats1.update(rate_stats2.mean, rate_stats2.pwr_sum_mean, rate_stats2.n, mask)

        assert numpy.allclose(rate_stats1.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats1.var, data.var(axis=0,))

    def test_with_mask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        n_vals = numpy.ones((nbins, nbins), numpy.uint)
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats = StreamingStats2D(nbins)
        for di, d in enumerate(data):
            rate_stats.update(d, d**2, n_vals, mask[di])

        data_masked = numpy.ma.array(data, mask=mask)
        rs_mean = numpy.ma.masked_where(numpy.isnan(rate_stats.mean), rate_stats.mean)
        rs_var = numpy.ma.masked_where(numpy.isnan(rate_stats.var), rate_stats.var)
        assert numpy.ma.allclose(rs_mean, data_masked.mean(axis=0))

        assert numpy.ma.allclose(rs_var, data_masked.var(axis=0,))

    def test_with_mask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins, nbins))
        n_vals = numpy.ones((nbins, nbins), numpy.uint)
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats1 = StreamingStats2D(nbins)
        for di, d in enumerate(data[:50]):
            rate_stats1.update(d, d**2, n_vals, mask[di])

        rate_stats2 = StreamingStats2D(nbins)
        for di, d in enumerate(data[50:]):
            rate_stats2.update(d, d**2, n_vals, mask[di])

        nomask = numpy.zeros((nbins,nbins), numpy.uint8)
        rate_stats1.update(rate_stats2.mean, rate_stats2.pwr_sum_mean, rate_stats2.n, nomask) 

        data_masked = numpy.ma.array(data, mask=mask)
        rs_mean = numpy.ma.masked_where(numpy.isnan(rate_stats1.mean), rate_stats1.mean)
        rs_var = numpy.ma.masked_where(numpy.isnan(rate_stats1.var), rate_stats1.var)

        assert numpy.ma.allclose(rs_mean, data_masked.mean(axis=0))
        assert numpy.ma.allclose(rs_var, data_masked.var(axis=0,))


class TestStreamingStats1D:
    def test_nomask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        n_vals = numpy.ones((nbins,), numpy.uint)
        mask = numpy.zeros((nbins,), numpy.uint8)

        rate_stats = StreamingStats1D(nbins)
        for d in data:
            rate_stats.update(d, d**2, n_vals, mask)

        assert numpy.allclose(rate_stats.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats.var, data.var(axis=0,))

    def test_nomask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        n_vals = numpy.ones((nbins,), numpy.uint)
        mask = numpy.zeros((nbins,), numpy.uint8)

        rate_stats1 = StreamingStats1D(nbins)
        for d in data[:50]:
            rate_stats1.update(d, d**2, n_vals, mask)

        rate_stats2 = StreamingStats1D(nbins)
        for d in data[50:]:
            rate_stats2.update(d, d**2, n_vals, mask)

        rate_stats1.update(rate_stats2.mean, rate_stats2.pwr_sum_mean, rate_stats2.n, mask)

        assert numpy.allclose(rate_stats1.mean, data.mean(axis=0))
        assert numpy.allclose(rate_stats1.var, data.var(axis=0,))

    def test_with_mask(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets, nbins))
        n_vals = numpy.ones((nbins,), numpy.uint)
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats = StreamingStats1D(nbins)
        for di, d in enumerate(data):
            rate_stats.update(d, d**2, n_vals, mask[di])

        data_masked = numpy.ma.array(data, mask=mask)
        rs_mean = numpy.ma.masked_where(numpy.isnan(rate_stats.mean), rate_stats.mean)
        rs_var = numpy.ma.masked_where(numpy.isnan(rate_stats.var), rate_stats.var)
        assert numpy.ma.allclose(rs_mean, data_masked.mean(axis=0))

        assert numpy.ma.allclose(rs_var, data_masked.var(axis=0,))

    def test_with_mask_groups(self):
        nbins = 100
        nsets = 10
        data = numpy.random.normal(size=(nsets,nbins))
        n_vals = numpy.ones((nbins,), numpy.uint)
        mask = numpy.random.randint(2, size=data.shape).astype(numpy.uint8)

        rate_stats1 = StreamingStats1D(nbins)
        for di, d in enumerate(data[:50]):
            rate_stats1.update(d, d**2, n_vals, mask[di])

        rate_stats2 = StreamingStats1D(nbins)
        for di, d in enumerate(data[50:]):
            rate_stats2.update(d, d**2, n_vals, mask[di])

        nomask = numpy.zeros((nbins,), numpy.uint8)
        rate_stats1.update(rate_stats2.mean, rate_stats2.pwr_sum_mean, rate_stats2.n, nomask) 

        data_masked = numpy.ma.array(data, mask=mask)
        rs_mean = numpy.ma.masked_where(numpy.isnan(rate_stats1.mean), rate_stats1.mean)
        rs_var = numpy.ma.masked_where(numpy.isnan(rate_stats1.var), rate_stats1.var)

        assert numpy.ma.allclose(rs_mean, data_masked.mean(axis=0))
        assert numpy.ma.allclose(rs_var, data_masked.var(axis=0,))
