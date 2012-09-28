from __future__ import division,print_function; __metaclass__ = type
import numpy
import west

from itertools import izip

from west.kinetics._kinetics import flux_assign, pop_assign, calc_rates #@UnresolvedImport

class RateAverager:
    '''Calculate bin-to-bin kinetic properties (fluxes, rates, populations) at
    1-tau resolution'''
    
    def __init__(self, bin_mapper, system=None, data_manager = None):
        self.bin_mapper = bin_mapper
        self.data_manager = data_manager or west.rc.get_data_manager()
        self.system = system or west.rc.get_system_driver()
        
    def calculate(self, iter_start=None, iter_stop=None):
        '''Read the HDF5 file and collect flux matrices and population vectors
        for each bin for each iteration in the range [iter_start, iter_stop)'''
        iter_start = iter_start or 1
        iter_stop = iter_stop or self.data_manager.current_iteration
        
        itercount = iter_stop - iter_start
        nbins = self.bin_mapper.nbins
        
        flux_matrices = numpy.zeros((itercount,nbins,nbins),numpy.float64)
        population_vectors = numpy.zeros((itercount,nbins), numpy.float64)

        pcoord_len = self.system.pcoord_len
        assign = self.bin_mapper.assign
        
        for iiter, n_iter in enumerate(xrange(iter_start, iter_stop)):
            iter_group = self.data_manager.get_iter_group(n_iter)
            
            # First, account for the flux due to recycling
            # We access the HDF5 file directly to avoid nearly 50% overhead of creating a ton of 
            # tiny NewWeightEntry objects            
            try:
                nwgroup = iter_group['new_weights']
            except KeyError:
                # No new weight data
                pass
            else:       
                index = nwgroup['index'][...]
                weights = index['weight']
                prev_init_pcoords = nwgroup['prev_init_pcoord'][...]
                new_init_pcoords = nwgroup['new_init_pcoord'][...]
                                    
                prev_init_assignments = assign(prev_init_pcoords)
                new_init_assignments  = assign(new_init_pcoords)
                
                flux_assign(weights, prev_init_assignments, new_init_assignments, flux_matrices[iiter])
                #for (weight,i,j) in izip (weights, prev_init_assignments, new_init_assignments):
                #    flux_matrices[iiter,i,j] += weight                
                del index                                
                del prev_init_pcoords, new_init_pcoords, prev_init_assignments, new_init_assignments, weights

                
            iter_group = self.data_manager.get_iter_group(n_iter)
            weights = iter_group['seg_index']['weight']
            
            initial_pcoords = iter_group['pcoord'][:,0]
            final_pcoords   = iter_group['pcoord'][:,pcoord_len-1]
            
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
        
        # Store references time series    
        self.flux_matrices = flux_matrices
        self.population_vectors = population_vectors
        
        # Calculate averages and standard deviations
        # These quantities are well-defined
        self.average_flux = flux_matrices.mean(axis=0)
        self.stdev_flux   = flux_matrices.std(axis=0)
        self.stderr_flux  = self.stdev_flux / len(self.flux_matrices)
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
        self.stdev_rate   = masked_rate_matrices.std(axis=0)
        self.stderr_rate  = self.stdev_rate / numpy.sqrt(numpy.sum(~masked_rate_matrices.mask,axis=0))
        
             
if __name__ == '__main__':
    # Tests this file on the west.h5 data in the current (sim root) directory    
    west.rc.read_config()
    system = west.rc.get_system_driver()
    data_manager = west.rc.get_data_manager()
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
    
    
    
        
        