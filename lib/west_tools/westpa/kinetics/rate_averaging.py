# Copyright (C) 2013 Matthew C. Zwier, Joshua L. Adelman, and Lillian T. Chong
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

import numpy
import westpa

from itertools import zip_longest
from collections import namedtuple

from westpa.kinetics._kinetics import flux_assign, pop_assign, calc_rates, StreamingStats1D, StreamingStats2D  #@UnresolvedImport


# Named tuple proxy for StreamingStats class
StreamingStatsTuple = namedtuple('StreamingStatsTuple', ['M1', 'M2', 'n'])


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def tuple2stats(stat_tuple):
    ndims = stat_tuple.M1.ndim
    assert ndims == 1 or ndims ==2

    if ndims == 2:
        stats = StreamingStats2D(stat_tuple.M1.shape)
    elif ndims == 1:
        stats = StreamingStats1D(stat_tuple.M1.shape[0])
    else:
        raise ValueError

    stats.M1 = stat_tuple.M1
    stats.M2 = stat_tuple.M2
    stats.n = stat_tuple.n

    return stats


def process_iter_chunk(bin_mapper, iter_indices, iter_data=None):
    '''Calculate the flux matrices and populations of a set of iterations specified
    by iter_indices. Optionally provide the necessary arrays to perform the calculation
    in iter_data. Otherwise get data from the data_manager directly.
    '''

    data_manager = westpa.rc.get_data_manager()
    system = westpa.rc.get_system_driver()

    itercount = len(iter_indices)
    nbins = bin_mapper.nbins

    flux_stats = StreamingStats2D((nbins, nbins))
    rate_stats = StreamingStats2D((nbins, nbins))
    pop_stats = StreamingStats1D(nbins)

    nomask1d = numpy.zeros((nbins,), numpy.uint8)
    nomask2d = numpy.zeros((nbins, nbins), numpy.uint8)
    rate_mask = numpy.zeros((nbins, nbins), numpy.uint8)

    flux_matrix = numpy.zeros((nbins, nbins), numpy.float64)
    rate_matrix = numpy.zeros((nbins, nbins), numpy.float64)
    population_vector = numpy.zeros((nbins,), numpy.float64)

    pcoord_len = system.pcoord_len
    assign = bin_mapper.assign

    for iiter, n_iter in enumerate(iter_indices):
        flux_matrix.fill(0.0)
        population_vector.fill(0.0)

        if iter_data:
            iter_group_name = 'iter_{:0{prec}d}'.format(int(n_iter), prec=data_manager.iter_prec)
            iter_group = iter_data[iter_group_name]
        else:
            iter_group = data_manager.get_iter_group(n_iter)

        # first, account for the flux due to recycling
        # we access the hdf5 file directly to avoid nearly 50% overhead of creating a ton of
        # tiny newweightentry objects
        try:
            nwgroup = iter_group['new_weights']
        except KeyError:
            # no new weight data
            pass
        else:
            if iter_data:
                index = None
                weights = nwgroup['weight']
                prev_init_pcoords = nwgroup['prev_init_pcoord']
                new_init_pcoords = nwgroup['new_init_pcoord']
            else:
                index = nwgroup['index'][...]
                weights = index['weight']
                prev_init_pcoords = nwgroup['prev_init_pcoord'][...]
                new_init_pcoords = nwgroup['new_init_pcoord'][...]

            prev_init_assignments = assign(prev_init_pcoords)
            new_init_assignments = assign(new_init_pcoords)

            flux_assign(weights, prev_init_assignments, new_init_assignments, flux_matrix)
            #for (weight,i,j) in izip (weights, prev_init_assignments, new_init_assignments):
            #    flux_matrices[iiter,i,j] += weight
            del index
            del prev_init_pcoords, new_init_pcoords, prev_init_assignments, new_init_assignments, weights

        #iter_group = data_manager.get_iter_group(n_iter)
        if iter_data:
            weights = iter_group['weight']
            initial_pcoords = iter_group['initial_pcoords']
            final_pcoords = iter_group['final_pcoords']
        else:
            weights = iter_group['seg_index']['weight']
            initial_pcoords = iter_group['pcoord'][:,0]
            final_pcoords = iter_group['pcoord'][:,pcoord_len-1]

        initial_assignments = assign(initial_pcoords)
        final_assignments = assign(final_pcoords)

        flux_assign(weights, initial_assignments, final_assignments, flux_matrix)
        pop_assign(weights, initial_assignments, population_vector)

        flux_stats.update(flux_matrix, nomask2d)
        pop_stats.update(population_vector, nomask1d)

        calc_rates(flux_matrix, population_vector, rate_matrix, rate_mask)
        rate_stats.update(rate_matrix, rate_mask)

        del weights
        del initial_assignments, final_assignments
        del initial_pcoords, final_pcoords
        del iter_group

    # Create namedtuple proxies for the cython StreamingStats objects
    # since the typed memoryviews class variables do not seem to return
    # cleanly from the zmq workers
    c_flux_stats = StreamingStatsTuple(flux_stats.M1, flux_stats.M2, flux_stats.n)
    c_rate_stats = StreamingStatsTuple(rate_stats.M1, rate_stats.M2, rate_stats.n)
    c_pop_stats = StreamingStatsTuple(pop_stats.M1, pop_stats.M2, pop_stats.n)

    return c_flux_stats, c_rate_stats, c_pop_stats 


class RateAverager():
    '''Calculate bin-to-bin kinetic properties (fluxes, rates, populations) at
    1-tau resolution'''

    def __init__(self, bin_mapper, system=None, data_manager=None, work_manager=None):
        self.bin_mapper = bin_mapper
        self.data_manager = data_manager or westpa.rc.get_data_manager()
        self.system = system or westpa.rc.get_system_driver()
        self.work_manager = work_manager or westpa.rc.get_work_manager()

    def extract_data(self, iter_indices):
        '''Extract data from the data_manger and place in dict mirroring the same
        underlying layout.'''

        data = {}
        pcoord_len = self.system.pcoord_len

        for n_iter in iter_indices:
            iter_group_name = 'iter_{:0{prec}d}'.format(int(n_iter), prec=self.data_manager.iter_prec)
            iter_group = self.data_manager.get_iter_group(n_iter)
            di = data[iter_group_name] = {}

            try:
                nwgroup = iter_group['new_weights']
            except KeyError:
                # no new weight data
                pass
            else:
                di_nw = di['new_weights'] = {}
                di_nw['weight'] = nwgroup['index'][...]['weight']
                di_nw['prev_init_pcoord'] = nwgroup['prev_init_pcoord'][...]
                di_nw['new_init_pcoord'] = nwgroup['new_init_pcoord'][...]

            di['weight'] = iter_group['seg_index']['weight']
            di['initial_pcoords'] = iter_group['pcoord'][:,0]
            di['final_pcoords'] = iter_group['pcoord'][:,pcoord_len-1]

        return data

    def task_generator(self, iter_start, iter_stop, block_size):
        for iter_block in grouper(block_size, range(iter_start, iter_stop)):
            iter_block = [x for x in iter_block if x is not None]
            iter_data = self.extract_data(iter_block)
            yield (process_iter_chunk, (self.bin_mapper, iter_block), {'iter_data': iter_data})
            del iter_data

    def calculate(self, iter_start=None, iter_stop=None, n_blocks=1, queue_size=1):
        '''Read the HDF5 file and collect flux matrices and population vectors
        for each bin for each iteration in the range [iter_start, iter_stop). Break
        the calculation into n_blocks blocks. If the calculation is broken up into
        more than one block, queue_size specifies the maxmimum number of tasks in
        the work queue.'''

        iter_start = iter_start or 1
        iter_stop = iter_stop or self.data_manager.current_iteration

        itercount = iter_stop - iter_start
        block_size = max(1, itercount // n_blocks)
        nbins = self.bin_mapper.nbins

        if n_blocks == 1:
            flux_stats_t, rate_stats_t, population_stats_t = process_iter_chunk(self.bin_mapper, list(range(iter_start, iter_stop)))

            flux_stats = tuple2stats(flux_stats_t)
            rate_stats = tuple2stats(rate_stats_t)
            population_stats = tuple2stats(population_stats_t)
        else:
            flux_stats = StreamingStats2D((nbins, nbins))
            rate_stats = StreamingStats2D((nbins, nbins))
            population_stats = StreamingStats1D(nbins)

            task_generator = self.task_generator(iter_start, iter_stop, block_size)

            for future in self.work_manager.submit_as_completed(task_generator, queue_size):
                chunk_flux_stats_t, chunk_rate_stats_t, chunk_pop_stats_t = future.get_result()

                chunk_flux_stats = tuple2stats(chunk_flux_stats_t)
                chunk_rate_stats = tuple2stats(chunk_rate_stats_t)
                chunk_pop_stats = tuple2stats(chunk_pop_stats_t)

                # Update statistics with chunked subsets
                flux_stats += chunk_flux_stats
                rate_stats += chunk_rate_stats
                population_stats += chunk_pop_stats

        self.average_flux = flux_stats.mean 
        self.stderr_flux = numpy.nan_to_num(numpy.sqrt(flux_stats.var) / flux_stats.n)

        self.average_populations = population_stats.mean 
        self.stderr_populations = numpy.nan_to_num(numpy.sqrt(population_stats.var) / population_stats.n)

        self.average_rate = rate_stats.mean 
        self.stderr_rate = numpy.nan_to_num(numpy.sqrt(rate_stats.var) / rate_stats.n)

        assert ~numpy.any(numpy.isinf(self.stderr_flux))
        assert ~numpy.any(numpy.isinf(self.stderr_rate))
        assert ~numpy.any(numpy.isinf(self.stderr_populations))


if __name__ == '__main__':
    # Tests this file on the west.h5 data in the current (sim root) directory
    westpa.rc.read_config()
    system = westpa.rc.get_system_driver()
    data_manager = westpa.rc.get_data_manager()
    data_manager.open_backing('r')
    averager = RateAverager(system.bin_mapper)
    averager.calculate()
    
    print('Population mean and standard error')
    print(averager.average_populations)
    print(averager.stderr_populations)
    
    print('Flux matrix, mean and standard error')
    print(averager.average_flux)
    print(averager.stderr_flux)
    
    print('Rate matrix, mean and standard error')
    print(averager.average_rate)
    print(averager.stderr_rate)

