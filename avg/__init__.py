'''
Set of tools for use with averaging together multiple WESTPA simulations.
The averaging code was written by Lewis Baker.
'''

import logging
log = logging.getLogger(__name__)

def average_data_slice(slice, dataset, iter, h5key):
    if len(slice.shape) = 2:
        return _average_2D_data_slice(slice, dataset, iter, h5key)
    if len(slice.shape) = 3:
        return _average_3D_data_slice(slice, dataset, iter, h5key)

def _average_2D_data_slice(slice, dataset, iter, h5key):
    # A function which takes in a slice, and a dataset, and adds on 
    # to it.
    slice[iter, :] += dataset[h5key][iter,:]
    return slice

def _average_3D_data_slice(slice, dataset, iter, h5key):
    # A function which takes in a slice, and a dataset, and adds on 
    # to it.
    slice[iter,:,:] += dataset[h5key][iter,:,:]
    return slice

def sigma_data_slice(slice, dataset, iter, h5key, ntrials, nstates, sigma = None):
    if len(slice.shape) = 2:
        return _create_sigma_2D_data(slice, dataset, iter, h5key, ntrials, nstates, sigma)
    if len(slice.shape) = 3:
        return _create_sigma_3D_data(slice, dataset, iter, h5key, ntrials, nstates, sigma)

def _create_sigma_3D_data(slice, dataset, iter, h5key, ntrials, nstates, sigma = None):
    if sigma = None:
        sigma = numpy.zeros( ( nstates, nstates ) )
    sigma = numpy.add( sigma, numpy.square( numpy.subtract( slice[iter,:,:],dataset[h5key][iter,:,:] ))/(ntrials - 1))
    return sigma

def _create_sigma_2D_data(slice, dataset, iter, h5key, ntrials, nstates, sigma = None):
    if sigma = None:
        sigma = numpy.zeros( ( nstates ) )
    sigma = numpy.add( sigma, numpy.square( numpy.subtract( slice[iter,:],dataset[h5key][iter,:] ))/(ntrials - 1))
    return sigma

def return_error_bounds(slice, sigma, i, niters, nstates):
    if len(slice.shape) = 2:
        return _return_error_bounds_3D_state(slice, sigma, i, niters, nstates)
    if len(slice.shape) = 3:
        return _return_error_bounds_3D_state(slice, sigma, i, niters, nstates)

def _return_error_bounds_3D_state(slice, sigma, i, j, niters, nstates):
    # Blank
    lbound = numpy.zeros( (niters, nstates, nstates) )
    ubound = numpy.zeros( (niters, nstates, nstates) )
    if sigma[i,j] == 0:
        lbound[iter,i,j] = 0
        ubound[iter,i,j] = 0
    else:
        bounds = scipy.stats.norm.interval(0.95, loc = slice[iter,i,j], scale = sigma[i,j] / numpy.sqrt(ntrials))
        c_radius = (bounds[1]-bounds[0])/2
        lbound[iter,i,j] = slice[iter,i,j] - c_radius
        ubound[iter,i,j] = slice[iter,i,j] + c_radius
    return ubound, lbound

def _return_error_bounds_2D_state(slice, sigma, i, niters, nstates):
    # Blank
    lbound = numpy.zeros( (niters, nstates) )
    ubound = numpy.zeros( (niters, nstates) )
    if sigma[i] == 0:
        lbound[iter,i] = 0
        ubound[iter,i] = 0
    else:
        bounds = scipy.stats.norm.interval(0.95, loc = slice[iter,i], scale = sigma[i] / numpy.sqrt(ntrials))
        c_radius = (bounds[1]-bounds[0])/2
        lbound[iter,i] = slice[iter,i] - c_radius
        ubound[iter,i] = slice[iter,i] + c_radius
    return ubound, lbound
