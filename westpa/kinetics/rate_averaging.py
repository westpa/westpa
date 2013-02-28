from __future__ import division,print_function; __metaclass__ = type
import numpy
import westpa

from itertools import izip, izip_longest

from westpa.kinetics._kinetics import flux_assign, pop_assign, calc_rates #@UnresolvedImport


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def process_iter_chunk(bin_mapper, iter_indices, iter_data=None):
    '''Calculate the flux matrices and populations of a set of iterations specified
    by iter_indices. Optionally provide the necessary arrays to perform the calculation
    in iter_data. Otherwise get data from the data_manager directly.
    '''

    data_manager = westpa.rc.get_data_manager()
    system = westpa.rc.get_system_driver()

    itercount = len(iter_indices)
    nbins = bin_mapper.nbins

    flux_matrices = numpy.zeros((itercount, nbins, nbins), numpy.float64)
    population_vectors = numpy.zeros((itercount, nbins), numpy.float64)

    pcoord_len = system.pcoord_len
    assign = bin_mapper.assign

    for iiter, n_iter in enumerate(iter_indices):
        if iter_data:
            iter_group_name = 'iter_{:0{prec}d}'.format(long(n_iter), prec=data_manager.iter_prec)
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

            flux_assign(weights, prev_init_assignments, new_init_assignments, flux_matrices[iiter])
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

        flux_assign(weights, initial_assignments, final_assignments, flux_matrices[iiter])
        pop_assign(weights, initial_assignments, population_vectors[iiter])
        #for (weight,i,j) in izip (weights, initial_assignments, final_assignments):
        #    flux_matrices[iiter,i,j] += weight
        #    population_vectors[iiter,i] += weight

        del weights
        del initial_assignments, final_assignments
        del initial_pcoords, final_pcoords
        del iter_group

    return flux_matrices, population_vectors, iter_indices


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
            iter_group_name = 'iter_{:0{prec}d}'.format(long(n_iter), prec=self.data_manager.iter_prec)
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
        for iter_block in grouper(block_size, xrange(iter_start, iter_stop)):
            iter_block = filter(lambda x: x is not None, iter_block)
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
            flux_matrices, population_vectors, _ = process_iter_chunk(self.bin_mapper, range(iter_start, iter_stop))
        else:
            flux_matrices = numpy.zeros((itercount,nbins,nbins),numpy.float64)
            population_vectors = numpy.zeros((itercount,nbins), numpy.float64)

            task_generator = self.task_generator(iter_start, iter_stop, block_size)

            for future in self.work_manager.submit_as_completed(task_generator, queue_size):
                chunk_flux_matrices, chunk_pop_vectors, iter_block = future.get_result()
                ii = numpy.array(iter_block) - iter_start
                flux_matrices[ii,:,:] = chunk_flux_matrices
                population_vectors[ii,:] = chunk_pop_vectors

        # store references time series
        self.flux_matrices = flux_matrices
        self.population_vectors = population_vectors

        # Calculate averages and standard deviations
        # These quantities are well-defined
        self.average_flux = flux_matrices.mean(axis=0)
        self.stdev_flux = flux_matrices.std(axis=0)
        self.stderr_flux = self.stdev_flux / len(self.flux_matrices)
        self.average_populations = population_vectors.mean(axis=0)
        self.stdev_populations = population_vectors.std(axis=0)
        self.stderr_populations = self.stdev_populations / len(self.flux_matrices)

        # Calculate rate matrices, while preparing masks for numpy masked arrays
        rate_matrices = numpy.zeros_like(flux_matrices)
        masks = numpy.ones(population_vectors.shape, numpy.bool_)
        calc_rates(flux_matrices, population_vectors, rate_matrices, masks)

        # Mask entries
        masked_rate_matrices = numpy.ma.array(rate_matrices)
        for rate_matrix, mask in izip(masked_rate_matrices, masks):
            rate_matrix[mask,mask] = numpy.ma.masked

        self.rate_matrices = masked_rate_matrices
        self.average_rate = masked_rate_matrices.mean(axis=0)
        self.stdev_rate = masked_rate_matrices.std(axis=0)
        self.stderr_rate = self.stdev_rate / numpy.sqrt(numpy.sum(~masked_rate_matrices.mask, axis=0))


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

