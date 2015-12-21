# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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

from __future__ import print_function, division; __metaclass__ = type
import logging

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csgraph
import h5py

from collections import Counter

import westpa
from west.data_manager import weight_dtype, n_iter_dtype
from westtools import (WESTMultiTool, WESTDataReader, IterRangeSelection,
                       ProgressIndicatorComponent)
from westpa import h5io
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

log = logging.getLogger('westtools.w_multi_reweight')

# TODO: remove unnecessary imports

def normalize(m):
    '''Row normalize the numpy array m, such that each row sums to one.  Return
    a new matrix.'''
    nm = m.copy()

    row_sum = m.sum(1)
    ii = np.nonzero(row_sum)[0]
    nm[ii,:] = m[ii,:] / row_sum[ii][:, np.newaxis]

    return nm

def steadystate_solve(K):
    '''Given a numpy array K representing a right stochastic matrix, use
    eigenvector/eigenvalue calculations to find the stationary distribution
    of the associated Markov chain in which probability is restricted to the
    largest, strongly-connected sub-graph. Return an array representing this
    probability vector.'''
    # Reformulate K to remove sink/source states
    n_components, component_assignments = csgraph.connected_components(K, connection="strong")
    largest_component = Counter(component_assignments).most_common(1)[0][0]
    components = np.where(component_assignments == largest_component)[0]

    ii = np.ix_(components, components)
    K_mod = K[ii]
    # This shouldn't really be needed any more, but shouldn't hurt.
    K_mod = normalize(K_mod)

    eigvals, eigvecs = np.linalg.eig(K_mod.T)
    eigvals = np.real(eigvals)
    eigvecs = np.real(eigvecs)

    maxi = np.argmax(eigvals)
    if not np.allclose(np.abs(eigvals[maxi]), 1.0):
        print('WARNING: Steady-state undetermined for current iteration')
        bin_prob = K.diagonal().copy()
        bin_prob = bin_prob / np.sum(bin_prob)
        return bin_prob

    sub_bin_prob = eigvecs[:, maxi] / np.sum(eigvecs[:, maxi])

    bin_prob = np.zeros(K.shape[0])
    bin_prob[components] = sub_bin_prob

    return bin_prob

def get_fluxmat_and_pop_sum(fluxmatH5_list, start_iter, stop_iter, nbins, 
                            recycling_bin_list=None):
    '''
    Find the sum of flux matrices and population vectors from a set of 
    simulations, for a given set of timepoints. 

    Returns: 
    (1) an N*N numpy array representing the sum of fluxes for each bin in all 
    simulations where the bin from which the flux originates is not a recycling
    bine. Elements for which no observations were made are set
    to ``NaN`` (not a number). 
    (2) an N*1 numpy array where each element represents the total probability 
    observed in each bin, for simulations where the bin is not a recycling bin. 

    Arguments:
    fluxmatH5_list: a list of flux matrix h5 files
    start_iter: integer; first weighted ensemble iteration to include in the
            average.
    stop_iter: integer; iteration at which to stop averaging.  Averages will 
            include UP TO stop_iter, BUT NOT stop_iter itself. 
    nbins: integer; the number of bins used in building the flux matrices. This
            should match the dimensions of the flux matrices in each HDF5 file 
            of fluxmatH5_list
    recycling_bin_list: A list of numpy arrays.  The order of the arrays should
            be the same as fluxmatH5_list. Each array should specify the indices
            of bins that were subject to recycling during the corresponding 
            weighted ensemble simulation. If all simulations were not subject
            to recycling (default), this option may be simply set to ``None``.
            If only select simulations were not subject to recycling, ``None`` 
            may be instead specified at the corresponding index of 
            recycling_bin_list. 
    '''
    # Get the population sums. First initialize the array to hold the sum
    pops_sum_arr = np.zeros(fluxmatH5_list[0]['bin_populations'].shape[1])
    temp_pop_arr = np.zeros(pops_sum_arr.shape)
    pops_sum_arr.fill(np.nan)
    temp_pop_arr.fill(np.nan)
    for isim, fluxmatH5 in enumerate(fluxmatH5_list):
        # Calculate the summed weight for this simulation.  Note that iteration 
        # indices are offset by 1 in 'bin_populations', being zero- rather than
        # one-indexed.
        temp_pop_arr = np.sum(fluxmatH5['bin_populations'][start_iter-1:stop_iter-1],
                              axis=0)
        # Set the recycling bins to a value of NaN (not a number)
        if recycling_bin_list is not None:
            if recycling_bin_list[isim] is not None:
                temp_pop_arr[recycling_bin_list[isim]] = np.nan 
        # Check the axis number here.
        pops_sum_arr = np.nansum(np.vstack((temp_pop_arr, pops_sum_arr)),
                                 axis=0)

    # Get the flux sums.  First initialize the array to hold the sum.
    flux_sum_arr = np.zeros((nbins,nbins))
    temp_flux_arr = np.zeros(flux_sum_arr.shape)
    flux_sum_arr.fill(np.nan)
    temp_flux_arr.fill(np.nan)
    obs_sum_arr = np.zeros(flux_sum_arr.shape)
    temp_obs_arr = np.zeros(flux_sum_arr.shape)

    for isim, fluxmatH5 in enumerate(fluxmatH5_list):
        rows = []
        cols = []
        obs = []
        flux = []
        # Nearly identical code to w_postanalysis_reweight
        for iiter in xrange(start_iter, stop_iter):
            iter_grp = fluxmatH5['iterations']['iter_{:08d}'.format(iiter)]
            # We need to reconstruct a dense matrix from the stored sparse
            # CooMatrix format.
            rows.append(iter_grp['rows'][...])
            cols.append(iter_grp['cols'][...])
            obs.append(iter_grp['obs'][...])
            flux.append(iter_grp['flux'][...])

        rows, cols, obs, flux = map(np.hstack, [rows, cols, obs, flux])
        # Convert to dense matrix AND sum the fluxes from each iteration, all 
        # in one line.
        temp_flux_arr = sp.coo_matrix((flux, (rows, cols)), shape=(nbins, nbins)).todense()
        # Convert the numpy "matrix" to a numpy "array", which can have more 
        # than 2 dimensions
        temp_flux_arr = np.array(temp_flux_arr)
        # Set the fluxes for the recycling bins (if any) to NaN
        if recycling_bin_list is not None:
            if recycling_bin_list[isim] is not None:
                temp_flux_arr[recycling_bin_list[isim]] = np.nan 
        # Add the fluxes for this simulation to the running total.
        flux_sum_arr = np.nansum(np.dstack((flux_sum_arr,temp_flux_arr)), axis=2)

        # Calculate the number of observations, accounting for the fact that we
        # ignore observations starting in recycling bins. Right now, this is not 
        # used anywhere. 
        temp_obs_arr = sp.coo_matrix((obs, (rows, cols)), shape=(nbins, nbins)).todense()
        temp_obs_arr[recycling_bin_list[isim]] = 0 
        obs_sum_arr += temp_obs_arr
    
    return flux_sum_arr, pops_sum_arr
        
    
def get_average_transition_mat_1(fluxmatH5_list, start_iter, stop_iter, nbins, 
                                 recycling_bin_list=None): 
    '''
    Find a transition matrix estimate from a set of simulations, for a given 
    set of timepoints. Return a numpy array representing the transition  
    matrix estimate, with elements for which no observations were made set
    to ``NaN`` (not a number). Estimate transition matrix elements as 
                          T_{i,j} = <w_{i,j}>/<w_i>
    where T_{i,j} represents the estimated probability of making a transition
    from bin i to bin j, w_{i,j} represents the observed flux from bin i to bin
    j in a given iteration of a given simulation, and w_i represents the
    probability observed in bin i during a given timepoint of a given 
    simulation. 

    Returns: 
    An N*N numpy array representing the estimated right-stochastic
    transition matrix.  Elements T_{i,j} of the matrix for which no probability
    was observed in bin i are set to "NaN" (not a number). 

    Arguments:
    fluxmatH5_list: a list of flux matrix h5 files
    start_iter: integer; first weighted ensemble iteration to include in the
            average.
    stop_iter: integer; iteration at which to stop averaging.  Averages will 
            include UP TO stop_iter, BUT NOT stop_iter itself. 
    nbins: integer; the number of bins used in building the flux matrices. This
            should match the dimensions of the flux matrices in each HDF5 file 
            of fluxmatH5_list
    recycling_bin_list: A list of numpy arrays.  The order of the arrays should
            be the same as fluxmatH5_list. Each array should specify the indices
            of bins that were subject to recycling during the corresponding 
            weighted ensemble simulation. If all simulations were not subject
            to recycling (default), this option may be simply set to ``None``.
            If only select simulations were not subject to recycling, ``None`` 
            may be instead specified at the corresponding index of 
            recycling_bin_list. 
    '''
    fluxmat_sum, pop_sum = get_fluxmat_and_pop_sum(
                                        fluxmatH5_list, start_iter, stop_iter, 
                                        nbins, 
                                        recycling_bin_list=recycling_bin_list
                                                    )
    
    # Finally, get the average transition matrix and return it.
    transmat = np.divide(fluxmat_sum, pop_sum[:,np.newaxis])
    return transmat


def cumulative_transmat_generator_1(fluxmatH5_list, start_iter, stop_iter, 
                                         step_iter, nbins, 
                                         recycling_bin_list=None):
    '''
    Generator for fast calculations of cumulatively averaged transition matrix 
    estimates. Uses the "ratio of the averages" for the bin-to-bin fluxes and
    populations.

    Find a transition matrix estimate from a set of simulations, for a given 
    set of timepoints. Return a numpy matrix representing the transition matrix
    estimate, with elements for which no observations were made set
    to ``NaN`` (Not a Number). Estimate transition matrix elements as 
                          T_{i,j} = <w_{i,j}>/<w_i>
    where T_{i,j} represents the estimated probability of making a transition
    from bin i to bin j, w_{i,j} represents the observed flux from bin i to bin
    j in a given iteration of a given simulation, and w_i represents the
    probability observed in bin i during a given timepoint of a given 
    simulation. 

    Yields: 
    N*N numpy arrays representing the estimated right-stochastic transition
    matrix based on data from windows of increasing width, starting with
    start_iter:start_iter+step_iter, and increasing in width by intervals of 
    step_iter until the window is start_iter:stop_iter. Elements T_{i,j} of 
    the matrix for which no probability was observed in bin i are set to "NaN" 
    (not a number). 

    Arguments:
    fluxmatH5_list: a list of flux matrix h5 files
    start_iter: integer; first weighted ensemble iteration to include in the
            average.
    stop_iter: integer; iteration at which to stop averaging.  Averages will 
            include UP TO stop_iter, BUT NOT stop_iter itself. 
    nbins: integer; the number of bins used in building the flux matrices. This
            should match the dimensions of the flux matrices in each HDF5 file 
            of fluxmatH5_list
    recycling_bin_list: A list of numpy arrays.  The order of the arrays should
            be the same as fluxmatH5_list. Each array should specify the indices
            of bins that were subject to recycling during the corresponding 
            weighted ensemble simulation. If all simulations were not subject
            to recycling (default), this option may be simply set to ``None``.
            If only select simulations were not subject to recycling, ``None`` 
            may be instead specified at the corresponding index of 
            recycling_bin_list. 
    '''
    # Build first cumulative_sum, from start_iter to start_iter + step_iter 
    # CHANGE: raise a more descriptive/accurate error here.
    if start_iter + step_iter > stop_iter:
        raise ValueError
    fluxmat_sum, pop_sum = get_fluxmat_and_pop_sum(fluxmatH5_list, start_iter,
                                           start_iter+step_iter, nbins, 
                                           recycling_bin_list=recycling_bin_list
                                                   )
    transmat = np.divide(fluxmat_sum, pop_sum[:,np.newaxis])
    # Yield the mean transition matrix from the first averaging window.
    yield transmat
    
    # From here on out, we can iterate.
    for cumsum_stop_iter in xrange(start_iter+step_iter, stop_iter, step_iter):
        # Get the sums between cumsum_stop_iter-step_iter and cumsum_stop_iter 
        temp_fluxmat_sum, temp_pop_sum = get_fluxmat_and_pop_sum(
                                         fluxmatH5_list, 
                                         cumsum_stop_iter-step_iter, 
                                         cumsum_stop_iter, nbins, 
                                         recycling_bin_list=recycling_bin_list
                                                                 )
        fluxmat_sum = np.nansum(np.dstack((fluxmat_sum, temp_fluxmat_sum)),
                                axis=2)
        print("Nonzero elements of fluxmat_sum: {:d}".format(np.count_nonzero(fluxmat_sum)))
        pop_sum += temp_pop_sum
        
        transmat = np.divide(fluxmat_sum, pop_sum[:,np.newaxis])

        # Yield the mean transition matrix from the  averaging window.
        yield transmat


def get_transmat_and_obsmat_sum(fluxmatH5_list, start_iter, stop_iter, nbins, 
                                recycling_bin_list=None):
    '''
    Find the sum of transition matrices from a set of simulations, for a given 
    set of timepoints. 

    Returns: 
    (1) a numpy array representing the sum of the transition matrices for each 
    simulation, for each timepoint for which the transition probability is 
    well-defined (ie, if the bin in which the transition originates is not 
    subject to recycling for the given simulation, and if this bin is not 
    unoccupied. 
    Elements for which no observations were made are set
    to ``NaN`` (not a number). 
    (2) a numpy array where each element represents the count of all timepoints
    for all simulations in which the transition probably represented in (1) was
    well-defined.

    Arguments:
    fluxmatH5_list: a list of flux matrix h5 files
    start_iter: integer; first weighted ensemble iteration to include in the
            average.
    stop_iter: integer; iteration at which to stop averaging.  Averages will 
            include UP TO stop_iter, BUT NOT stop_iter itself. 
    nbins: integer; the number of bins used in building the flux matrices. This
            should match the dimensions of the flux matrices in each HDF5 file 
            of fluxmatH5_list
    recycling_bin_list: A list of numpy arrays.  The order of the arrays should
            be the same as fluxmatH5_list. Each array should specify the indices
            of bins that were subject to recycling during the corresponding 
            weighted ensemble simulation. If all simulations were not subject
            to recycling (default), this option may be simply set to ``None``.
            If only select simulations were not subject to recycling, ``None`` 
            may be instead specified at the corresponding index of 
            recycling_bin_list. 
    '''
    
    # Array for holding the running sum of transition matrix
    transmat_sum = np.empty((nbins, nbins), np.float64)
    transmat_sum.fill(np.nan)
    temp_array = np.zeros((nbins, nbins), np.float64)
    # Array for holding the count of elements added to the transmat_sum 
    obsmat_sum = np.zeros((nbins, nbins), np.int64)
    
    # Load all the populations into physical memory ahead of time.  This should
    # be small enough to store without running into memory issues. Only read in
    # the specified iteration range. Note that the iteration indexing of these 
    # vectors is offset from the h5file itself!
    pops_list = []
    for fluxmatH5 in fluxmatH5_list:
        pops_list.append(
                np.array(fluxmatH5['bin_populations'][start_iter-1:stop_iter-1, :])
                         ) 

    np.seterr(divide='ignore')
    # Iterate over each iteration of each simulation
    for iiter in xrange(start_iter, stop_iter):
        for isim, fluxmatH5 in enumerate(fluxmatH5_list):
            iter_group = fluxmatH5['iterations/iter_{:08d}'.format(iiter)]
            # Reset the temporary matrix to all zeros.
            temp_array[...] = 0 
            # Load flux matrix from the HDF5 file. The matrix is stored based
            # on the ``coo_matrix`` format from the sparse module of the SciPy 
            # library. The row index, column index, and associated data is 
            # stored in vectors, for nonzero elements of the array.
            row_idxs = iter_group['rows']
            col_idxs = iter_group['cols']
            flux_data = iter_group['flux']
            # Reconstruct the (dense) matrix.
            temp_array[row_idxs, col_idxs] = flux_data
            # Divide each row of the matrix by the population in corresponding
            # bin.  This gives ``NaN`` for values where the population of the 
            # corresponding bin is zero (ie, 0/0). 
            temp_array = np.divide(temp_array, 
                                   pops_list[isim][iiter-start_iter][:,np.newaxis]
                                   )
            # Set rows corresponding to bins involved in recycling to NaN, as 
            # well. 
            if recycling_bin_list is not None:
                if recycling_bin_list[isim] is not None:
                    temp_array[recycling_bin_list[isim]] = np.nan 
            # Make a mask that lets through the values that are not 
            # ``Not a Number``
            good_value_mask = np.isfinite(temp_array) 
            # Add the good values to the running sum
            # We could also experiment with weighting elements by the number of
            # observations here!
            transmat_sum = np.nansum(np.dstack((transmat_sum, temp_array)),
                                     axis=2)
            obsmat_sum += good_value_mask
    return transmat_sum, obsmat_sum

def get_average_transition_mat_2(fluxmatH5_list, start_iter, stop_iter, nbins, 
                                 recycling_bin_list=None): 
    '''
    Find an average transition matrix from a set of simulations, for a given 
    set of timepoints. Return a numpy array representing the average 
    transition matrix, with elements for which no observations were made set
    to ``NaN`` (not a number). 

    Returns: 
    An N*N numpy array representing the estimated right-stochastic
    transition matrix.  Elements T_{i,j} of the matrix for which no probability
    was observed in bin i are set to "NaN" (not a number). 

    Arguments:
    fluxmatH5_list: a list of flux matrix h5 files
    start_iter: integer; first weighted ensemble iteration to include in the
            average.
    stop_iter: integer; iteration at which to stop averaging.  Averages will 
            include UP TO stop_iter, BUT NOT stop_iter itself. 
    nbins: integer; the number of bins used in building the flux matrices. This
            should match the dimensions of the flux matrices in each HDF5 file 
            of fluxmatH5_list
    recycling_bin_list: A list of numpy arrays.  The order of the arrays should
            be the same as fluxmatH5_list. Each array should specify the indices
            of bins that were subject to recycling during the corresponding 
            weighted ensemble simulation. If all simulations were not subject
            to recycling (default), this option may be simply set to ``None``.
            If only select simulations were not subject to recycling, ``None`` 
            may be instead specified at the corresponding index of 
            recycling_bin_list. 
    '''
    transmat_sum, obsmat_sum = get_transmat_and_obsmat_sum(
                                        fluxmatH5_list, start_iter, stop_iter, 
                                        nbins, 
                                        recycling_bin_list=recycling_bin_list
                                                           )
    
    # Finally, get the average transition matrix and return it.
    transmat = np.divide(transmat_sum, obsmat_sum)
    return transmat
            
def cumulative_transmat_generator_2(fluxmatH5_list, start_iter, stop_iter, 
                                         step_iter, nbins, 
                                         recycling_bin_list=None):
    '''
    Generator for fast calculations of cumulatively averaged transition matrix 
    estimates. Uses the "average ratio" of the bin-to-bin fluxes and the bin 
    populations.

    Find an average transition matrix from a set of simulations, for a given 
    set of timepoints. Return a numpy matrix representing the average 
    transition matrix, with elements for which no observations were made set
    to ``NaN`` (Not a Number). 

    Yields: 
    N*N numpy arrays representing the estimated right-stochastic transition
    matrix based on data from windows of increasing width, starting with
    start_iter:start_iter+step_iter, and increasing in width by intervals of 
    step_iter until the window is start_iter:stop_iter. Elements T_{i,j} of 
    the matrix for which no probability was observed in bin i are set to "NaN" 
    (not a number). 

    Arguments:
    fluxmatH5_list: a list of flux matrix h5 files
    start_iter: integer; first weighted ensemble iteration to include in the
            average.
    stop_iter: integer; iteration at which to stop averaging.  Averages will 
            include UP TO stop_iter, BUT NOT stop_iter itself. 
    step_iter: for each cumulative mean, include another step_iter iterations
    nbins: integer; the number of bins used in building the flux matrices. This
            should match the dimensions of the flux matrices in each HDF5 file 
            of fluxmatH5_list
    recycling_bin_list: A list of numpy arrays.  The order of the arrays should
            be the same as fluxmatH5_list. Each array should specify the indices
            of bins that were subject to recycling during the corresponding 
            weighted ensemble simulation. If all simulations were not subject
            to recycling (default), this option may be simply set to ``None``.
            If only select simulations were not subject to recycling, ``None`` 
            may be instead specified at the corresponding index of 
            recycling_bin_list. 
    '''
    # Build first cumulative_sum, from start_iter to start_iter + step_iter 
    # CHANGE: raise a more descriptive/accurate error here.
    if start_iter + step_iter >= stop_iter:
        raise ValueError
    transmat_sum, obsmat_sum = get_transmat_and_obsmat_sum(
                                        fluxmatH5_list, start_iter, 
                                        start_iter+step_iter, nbins, 
                                        recycling_bin_list=recycling_bin_list
                                                           )
    transmat = np.divide(transmat_sum, obsmat_sum)
    # Yield the mean transition matrix from the first averaging window.
    yield transmat
    
    # From here on out, we can iterate.
    for cumsum_stop_iter in xrange(start_iter+step_iter, stop_iter, step_iter):
        # Get the sums between cumsum_stop_iter-step_iter and cumsum_stop_iter 
        temp_transmat_sum, temp_obsmat_sum = get_transmat_and_obsmat_sum(
                                        fluxmatH5_list, 
                                        cumsum_stop_iter-step_iter, 
                                        cumsum_stop_iter, nbins, 
                                        recycling_bin_list=recycling_bin_list
                                                                         )
        transmat_sum = np.nansum(np.dstack((transmat_sum, temp_transmat_sum)),
                                 axis=2)
        obsmat_sum += temp_obsmat_sum

        transmat = np.divide(transmat_sum, obsmat_sum)

        # Yield the mean transition matrix from the  averaging window.
        yield transmat


def reweight(transition_matrix, nstates, nbins, state_map): 
    ''' Apply the non-Markovian reweighting procedure to the given transition 
    matrix.'''
    # TODO: add in observation threshold?

    # Need to get rid of "NaN" values in the trans_mat.  Set these to zero.
    transition_matrix[np.isnan(transition_matrix)] = 0.0
    # Solve for the steady-state (stationary) distribution.
    rw_bin_probs = steadystate_solve(transition_matrix)

    bin_last_state_map = np.tile(np.arange(nstates, dtype=np.int), nbins)
    bin_state_map = np.repeat(state_map[:-1], nstates)

    rw_color_probs = np.bincount(bin_last_state_map, weights=rw_bin_probs) 
    rw_state_probs = np.bincount(bin_state_map, weights=rw_bin_probs)

    rw_bin_transition_matrix = transition_matrix

    ii = np.nonzero(transition_matrix)

    rw_state_flux = calc_state_flux(rw_bin_transition_matrix[ii], ii[0], ii[1], rw_bin_probs, 
            bin_last_state_map, bin_state_map, nstates)

    return rw_state_flux, rw_color_probs, rw_state_probs, rw_bin_probs, rw_bin_transition_matrix


def calc_state_flux(trans_matrix, index1, index2, bin_probs, bin_last_state_map, bin_state_map, nstates):
    state_flux = np.zeros((nstates, nstates), np.float64)
    
    n_trans = index1.shape[0]
    for k in xrange(n_trans):
        ii = bin_last_state_map[index1[k]]
        jj = bin_state_map[index2[k]]

        if jj != nstates:
            state_flux[ii, jj] += trans_matrix[k] * bin_probs[index1[k]]

    return state_flux




class WMultiReweightTool(WESTMultiTool):
    prog ='w_multi_reweight'
    description = '''\
Calculate average rates from weighted ensemble data using the postanalysis
reweighting scheme. Bin assignments (usually "assignments.h5") and pre-calculated 
iteration flux matrices (usually "flux_matrices.h5") data files must have been 
previously generated using w_postanalysis_matrix.py (see "w_assign --help" and 
"w_kinetics --help" for information on generating these files).


-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, usually "kinrw.h5") contains the following
dataset:

  /state_prob_evolution [window,state]
    The reweighted state populations based on windows

  /color_prob_evolution [window,state]
    The reweighted populations last assigned to each state based on windows

  /bin_prob_evolution [window, bin]
    The reweighted populations of each bin based on windows. Bins contain
    one color each, so to recover the original un-colored spatial bins,
    one must sum over all states.

  /conditional_flux_evolution [window,state,state]
    (Structured -- see below). State-to-state fluxes based on windows of
    varying width
    
The structure of the final dataset is as follows:

  iter_start
    (Integer) Iteration at which the averaging window begins (inclusive).
    
  iter_stop
    (Integer) Iteration at which the averaging window ends (exclusive).
    
  expected
    (Floating-point) Expected (mean) value of the rate as evaluated within
    this window, in units of inverse tau.

(Optional)
If specified, transition matrices are saved (-t/--save-transition-matrices)
in /iterations/iter_08d/ as a sparse matrix.


-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WMultiReweightTool, self).__init__()
        
        self.data_reader = WESTDataReader()
        self.progress = ProgressIndicatorComponent()
        
        self.output_filename = None
        self.kinetics_filename = None
        self.assignment_filename = None
        
        self.output_file = None
        self.assignments_file = None
        self.kinetics_file = None
        
        self.evolution_mode = None
        
    def add_args(self, parser):
        self.progress.add_args(parser)
        self.data_reader.add_args(parser)
        ##self.iter_range.include_args['iter_step'] = True
        ##self.iter_range.add_args(parser)

        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('-y', '--yaml', dest='yamlpath', 
                             metavar='YAMLFILE', 
                             default='multi_reweight.yaml',
                             required=True,
                             help='''Load options from YAMLFILE. For each 
                             simulation, specify an assignments file and flux
                             matrices file. Search for files in 
                             ['simulations'][SIMNAME]['assignments'] and
                             ['simulations'][SIMNAME]['kinetics']. Also search
                             for lists of bins indicating those subject to 
                             recycling during the simulation 
                             ('recycling_bin_list'), or lists of states
                             indicating those subject to recycling during the 
                             simulation ('recycling_state_list').  Search in
                             ['simulations'][SIMNAME]['recycling_bin_list']
                             and ['simulations'][SIMNAME]['recycling_state_list'].
                             If supplied, elements of the transition rate matrix
                             corresponding to rates originating from bins in 
                             either of these lists are excluded from
                             contributing to the estimate for that element of
                             the transition rate matrix. 
                             ''')

        iogroup.add_argument('-o', '--output', dest='output', default='kinrw.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

        iogroup.add_argument('-t', '--save-transition-matrices', 
                             dest='save_transition_matrices', 
                             action='store_true',
                             help='''Include transition matrices in output file. 
                             Save transition matrices in ['iterations/iter_08d/]
                             as three datasets keyed as 'rows', 'cols', and 'k',
                             corresponding to a sparse matrix format. The  
                             transition matrix ``T`` may be reconstructed as 
                             T[rows][cols] = k, with T elsewhere zero.
                             Transition matrices are right stochastic, with the
                             i'th, j'th element corresponding to the transition
                             rate constant estimate from bin i to bin j for the
                             supplied lag time.''')

        cogroup = parser.add_argument_group('calculation options')
        cogroup.add_argument('-e', '--evolution-mode', choices=['cumulative', 'blocked'], default='cumulative',
                             help='''How to calculate time evolution of rate estimates.
                             ``cumulative`` evaluates rates over windows starting with --start-iter and getting progressively
                             wider to --stop-iter by steps of --step-iter.
                             ``blocked`` evaluates rates over windows of width --step-iter, the first of which begins at
                             --start-iter.''')
        cogroup.add_argument('--window-frac', type=float, default=1.0,
                             help='''Fraction of iterations to use in each window when running in ``cumulative`` mode.
                             The (1 - frac) fraction of iterations will be discarded from the start of each window.''')
        # Currently not implemented! (TODO) 
        #cogroup.add_argument('--obs-threshold', type=int, default=1,
        #                     help='''The minimum number of observed transitions between two states i and j necessary to include
        #                     fluxes in the reweighting estimate''')
        cogroup.add_argument('--first-iter', default='1', dest='first_iter',
                             type=int,
                             help='''Begin calculations starting at iteration 
                             FIRST_ITER''')
        cogroup.add_argument('--last-iter', default=None, dest='last_iter',
                             type=int,
                             help='''Include iterations up to and including 
                             LAST_ITER in calculations.''')
        cogroup.add_argument('--step-iter', default=None, dest='step_iter',
                             type=int,
                             help='''Move averaging windows forward in widths
                             of STEP_ITER.  For cumulative average, begin
                             averaging with the first window including
                             iterations from FIRST_ITER to FIRST_ITER + 
                             STEP_ITER, and then expand with window forward in 
                             widths of STEP_ITER.  For blocked averaging, use 
                             blocks of size STEP_ITER; the first window will 
                             include data from iteration FIRST_ITER to 
                             FIRST_ITER + STEP_ITER, while the second iteration 
                             will include data from iteration FIRST_ITER + 
                             STEP_ITER to FIRST_ITER + 2*STEP_ITER, etc.''')
        cogroup.add_argument('-s', '--estimation-scheme', choices=[1,2],
                             dest='estimation_scheme',
                             default=1, type=int,
                             help='''Estimate bin-to-bin transition rates using
                             one of the following schemes. 

                             Scheme (1) [default]: use the ratio of the average 
                             flux and average bin population, i.e. 
                             <w_{i,j}>/<w_i>, where w_{i,j} is the flux from bin
                             i to bin j in a given iteration, and w_i is the 
                             probability in bin i in a given iteration.

                             Scheme (2): use the average ratio of flux and bin
                             population, i.e. <w_{i,j}/w_i>.''') 
                              
                                                         
    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        #with self.data_reader:
        ##self.iter_range.process_args(args, default_iter_step=None)
        #if self.iter_range.iter_step is None:
        self.iter_range = {'first_iter': args.first_iter,
                           'last_iter':  args.last_iter,
                           'step_iter':  args.step_iter  
                           }
    
        
        self.output_filename = args.output
        self.evolution_mode = args.evolution_mode
        self.estimation_scheme = args.estimation_scheme
        self.evol_window_frac = args.window_frac
        if self.evol_window_frac <= 0 or self.evol_window_frac > 1:
            raise ValueError('Parameter error -- fractional window defined by '
                             '--window-frac must be in (0,1]')
        if args.save_transition_matrices:
            self.save_transition_matrices = True
        else:
            self.save_transition_matrices = False

        # Load the yaml input file; Make self.yamlargdict available
        self.parse_from_yaml(args.yamlpath)
        # Open indicated files and check for consistency.
        self.open_files()
        self.check_consistency_of_input_files()
        self.load_recycling_bins()


    def load_recycling_bins(self):
        '''Load the list of bins that were subject to recycling in each
        simulation.  Check in ``yamldict`` under ['simulations'][SIMNAME],
        looking for the keys 'recycling_bin_list' or 'recycling_state_list'.
        Use the state_map in each assignments file to convert state indices to 
        bin indices, if need be, and union over the bin indices given in 
        'recycling_bin_list' and those corresponding to states indicated in 
        'recycling_state_list'.'''
        # Get the list of recycling bins from the YAML file, setting to None
        # if not specified.
        self.recycling_bin_list = [] 
        for isim, simname in enumerate(self.yamlargdict['simulations'].keys()):
            simdict = self.yamlargdict['simulations'][simname]
            binset = set()
            if 'recycling_bin_list' in simdict.keys(): 
                binset = binset.union(set(simdict['recycling_bin_list'])) 
            if 'recycling_state_list' in simdict.keys():
                state_idxs = simdict['recycling_state_list']
                state_map = self.assignments_file_list[isim]['state_map']
                bin_idxs = set() 
                for state_idx in state_idxs:
                    bin_idxs = bin_idxs.union( 
                                   set(np.where(state_map[...] == state_idx)[0])
                                              )
                binset = binset.union(bin_idxs)
            if binset:
                self.recycling_bin_list.append(list(binset))
            else:
                self.recycling_bin_list.append(None)


    def adjust_recycling_bins_for_color(self, nstates):
        '''For all the specified recycling bins, adjust the indices to match
        the colored matrix (with int ``nstates`` colors) indices.  E.g., if 
        the user specifies bin index ``i`` as a recycling bin, and there are
        ``nstates`` different colors, then indices 
            nstates*i
            nstates*i + 1
            nstates*i + 2 
            ...
            nstates*i + (nstates-1)
        correspond to recycling bins in the colored matrix.'''
        for isim, recycling_bins in enumerate(self.recycling_bin_list):
            if recycling_bins is not None:
                new_recycling_bins = []
                for bin_index in recycling_bins:
                    for j in xrange(nstates):
                       new_recycling_bins.append(nstates*bin_index + j)
                self.recycling_bin_list[isim] = new_recycling_bins 
            
         
    def open_files(self):
        '''
        Create the output file.  Open input files, including assignments and 
        kinetics files for each simulation specified.  Make the output file 
        available as an attributue of ``self``, ``output_file``.  Similarly,
        make lists of the input files available as attributes, 
        ``assignments_file_list`` and ``kinetics_files_list``.
        '''
        # What is opening this file before I do?
        # Need to figure out what's going on here...
        self.output_file = h5io.WESTPAH5File(self.output_filename, 'a', 
                                             creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.assignments_file_list = []
        self.kinetics_file_list = []
        self.simname_list = []
        # Open the assignments and kinetics (flux matrices) files for each
        # simulation.  assignments_file_list[i] and kinetcs_file_list[i] 
        # correspond to the same simulation.
        for simname in self.yamlargdict['simulations'].keys():
            assignments_filename = self.yamlargdict['simulations'][simname]\
                                                   ['assignments']
            assignH5 = h5io.WESTPAH5File(assignments_filename, 'r')
            # Double check that adding this file to the list adds the file
            # itself, and not a pointer to the assignH5 variable.
            self.assignments_file_list.append(assignH5)  
            kinetics_filename = self.yamlargdict['simulations'][simname]\
                                                ['kinetics']
            kinetH5 = h5io.WESTPAH5File(kinetics_filename, 'r')
            self.kinetics_file_list.append(kinetH5)

            # Keep track of the name associated with each simulation. This tool
            # uses these names later to help give descriptive errors!
            self.simname_list.append(simname)
            # Check that the specified assignment and kinetics (flux matrices)
            # files span the requested iterations.

    def check_consistency_of_input_files(self):
        '''
        Check the input files for consistency.  This requires that the number of
        bins used to create all the assignment files and all the kinetics (flux
        matrices) files is the same. Additionally, "npts" (timepoints per
        iteration) must be consistent between flux matrix files.
        '''
        old_nrows = None
        for i, fmH5 in enumerate(self.kinetics_file_list):
            nrows = fmH5.attrs['nrows'] 
            if old_nrows is not None:
                if nrows != old_nrows:
                    prev_simname = self.simname_list[i-1]
                    curr_simname = self.simname_list[i]
                    # Could make this error more descriptive than ValueError
                    raise ValueError("Number of bins specified in simulation "
                                     "{:s} ({:d}) is not the same as simulation"
                                     " {:s} ({:d}).".format(prev_simname, 
                                                            old_nrows,
                                                            curr_simname,
                                                            nrows)
                                     )
            old_nrows = nrows
        # Check that the number of bins used in making the assignments files
        # does not exceed the number of rows in the transition matrix; this 
        # could happen if the user specifies a kinetics file that was not 
        # created from the given assignments file.
        for i, assignH5 in enumerate(self.assignments_file_list):
            statecount = assignH5.attrs['nstates']
            bincount = assignH5.attrs['nbins']
            # Check that the number of bins used in the kinetics and assignment
            # files are consistent; valid only for non-markovian (colored) matrix
            if bincount*statecount != nrows:
                curr_simname = self.simname_list[i]
                raise ValueError("Number of bins ({:d}) used in assignments "
                                 "file for simulation {:s} does not match "
                                 "the number of rows used in kinetics "
                                 "files! ({:d})".format(bincount,
                                                        curr_simname,
                                                        nrows)
                                 )

        # Check that npts (timepoints per iteration) is consistent.
        prev_npts = None
        for i, fmH5 in enumerate(self.kinetics_file_list):
            npts = fmH5.attrs['npts']
            if prev_npts is not None:
                if npts != prev_npts:
                    curr_simname = self.simname_list[i]
                    prev_simname = self.simname_list[i-1]
                    raise ValueError("Number of timepoints ({:d}) used in "
                                     "kinetics (flux matrices) file for "
                                     "simulation {:s} does not match the number"
                                     " of timepoints used in simulation {:s} "
                                     "({:d})".format(npts, curr_simname,
                                                     prev_npts, prev_simname)
                                     ) 
            prev_npts = npts


    def get_iteration_ranges(self):
        '''Set default iteration ranges for those that were not specified.'''
        # Set the last iteration to be used in analysis to the least of the 
        # last iterations specified in any of the assignments or flux matrix 
        # files
        if self.iter_range['last_iter'] is None:
            least_last_iter = self.assignments_file_list[0]\
                    ['assignments'].attrs['iter_stop']-1 
            for isim, simname in enumerate(self.simname_list):
                # Check the "-1" term!
                curr_last_iter = self.assignments_file_list[isim]\
                        ['assignments'].attrs['iter_stop']-1
                if least_last_iter > curr_last_iter:
                    least_last_iter = curr_last_iter 
                curr_last_iter = self.kinetics_file_list[isim]\
                        ['bin_populations'].attrs['iter_stop']-1
                if least_last_iter > curr_last_iter:
                    least_last_iter = curr_last_iter 
            self.iter_range['last_iter'] = least_last_iter 
        # Use about 10 blocks by default, if iter_step is not specified.
        if self.iter_range['step_iter'] is None:
            self.iter_range['step_iter'] = max(1, 
              (self.iter_range['last_iter']-self.iter_range['first_iter']) // 10
                                               )


    def check_iteration_ranges(self):
        '''Check that the iteration ranges specified or calculated from the
        input files are actually available from the input files, and make 
        sense.'''
        if self.iter_range['last_iter'] <= self.iter_range['first_iter']:
            raise ValueError('Error. The specified or calculated iteration does'
                             ' not make sense (start_iter: {:d} >= '
                             'last_iter: {:d}). Check the iteration range you ' 
                             'specified (if applicable), or otherwise check the'
                             ' supplied input files to make sure they have '
                             'overlapping iteration ranges.')
        # Check that all the supplied files contain the needed iteration ranges
        for isim in xrange(len(self.simname_list)):
            assignH5 = self.assignments_file_list[isim]
            kinetH5 = self.kinetics_file_list[isim]
            if not (assignH5['assignments'].attrs['iter_stop'] >= \
                    self.iter_range['last_iter'] and \
                    assignH5['assignments'].attrs['iter_start'] <= \
                    self.iter_range['first_iter']):
                raise ValueError('Assignments data {:s} do not span the '
                                 'requested iterations ({:d} to {:d}) for '
                                 'simulation {:s}'.format(
                                         assignments_filename, 
                                         self.iter_range['first_iter'],
                                         self.iter_range['last_iter'], simname
                                                          )
                                 )

            if not (kinetH5['bin_populations'].attrs['iter_stop'] >= \
                    self.iter_range['last_iter'] and \
                    kinetH5['bin_populations'].attrs['iter_start'] <= \
                    self.iter_range['first_iter']):
                raise ValueError('Kinetics data {:s} do not span the '
                                 'requested iterations ({:d} to {:d}) for '
                                 'simulation {:s}'.format(
                                         assignments_filename, 
                                         self.iter_range['first_iter'],
                                         self.iter_range['last_iter'], simname
                                                          )
                                 )
             
    def go(self):
        '''Run the main analysis.'''
        # First, check over the input files and set up the iteration ranges. 
        self.get_iteration_ranges()
        self.check_iteration_ranges()
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            # Get the number of states and the number of bins for the
            # simulations (this should be the same across all the simulations,
            # since this was already checked with 
            # ``check_consistency_of_input_files``).
            nstates = self.assignments_file_list[0].attrs['nstates']
            nbins = self.assignments_file_list[0].attrs['nbins']

            # Adjust the recycling bins for the colored matrix
            self.adjust_recycling_bins_for_color(nstates)

            # Get the state labels for each simulation.  These should also be
            # consistent across simulations, but this is up to the user to make
            # sure of.  We could add in a check here later if we decide it
            # makes more sense.
            state_labels = self.assignments_file_list[0]['state_labels'][...]
            # State map *should* also be consistent, but we do not explicitly 
            # enforce that it is.
            state_map = self.assignments_file_list[0]['state_map'][...]
            # Get nfbins and npts (forced to be consistent by
            # ``check_consistency_of_input_files``).
            nfbins = self.kinetics_file_list[0].attrs['nrows']
            npts = self.kinetics_file_list[0].attrs['npts']

            assert nstates == len(state_labels)
            assert nfbins == nbins * nstates

            start_iter, stop_iter, step_iter = self.iter_range['first_iter'], \
                    self.iter_range['last_iter'], self.iter_range['step_iter']

            start_pts = range(start_iter, stop_iter, step_iter)
            flux_evol = np.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
            color_prob_evol = np.zeros((len(start_pts), nstates))
            state_prob_evol = np.zeros((len(start_pts), nstates))
            bin_prob_evol = np.zeros((len(start_pts), nfbins))
            pi.new_operation('Calculating flux evolution', len(start_pts))

            if self.save_transition_matrices:
                self.output_file.create_group('iterations')

            # Set up the generator, if needed; The generator reduces repetitive
            # calculations 
            if self.evolution_mode == 'cumulative' \
                    and self.evol_window_frac == 1.0:
                if self.estimation_scheme == 1:
                    gen = cumulative_transmat_generator_1(self.kinetics_file_list,
                                                          start_iter, stop_iter,
                                                          step_iter, nfbins,
                                                          self.recycling_bin_list
                                                          ) 
                if self.estimation_scheme == 2:
                    gen = cumulative_transmat_generator_2(self.kinetics_file_list,
                                                          start_iter, stop_iter,
                                                          step_iter, nfbins,
                                                          self.recycling_bin_list
                                                          ) 
            else:
                # set "get_transmat" to be the correct function for estimating
                # the transition matrix, based on the specified scheme.
                if self.estimation_scheme == 1:
                    get_transmat = get_average_transition_mat_1
                if self.estimation_scheme == 2:
                    get_transmat = get_average_transition_mat_2

            for iblock, start in enumerate(start_pts):
                pi.progress += 1

                stop = min(start + step_iter, stop_iter)
                # Get the transition matrix; depends on the averaging scheme 
                if self.evolution_mode == 'cumulative' \
                        and self.evol_window_frac == 1.0:
                    transition_matrix = gen.next() 
                else: 
                    if self.evolution_mode == 'cumulative':
                        windowsize = max(1, int(self.evol_window_frac * (stop - start_iter)))
                        block_start = max(start_iter, stop - windowsize)
                    else:   # self.evolution_mode == 'blocked'
                        block_start = start
                    transition_matrix = get_transmat(self.kinetics_file_list,
                                                     block_start, stop, 
                                                     nfbins,
                                                     self.recycling_bin_list)
                # Apply reweighting procedure 
                rw_state_flux, rw_color_probs, rw_state_probs, \
                        rw_bin_probs, rw_bin_flux = reweight(transition_matrix,
                                                             nstates, nbins,
                                                             state_map)
                for k in xrange(nstates):
                    for j in xrange(nstates):
                        # Normalize such that we report the flux per tau (tau being the weighted ensemble iteration)
                        # npts always includes a 0th time point
                        flux_evol[iblock]['expected'][k,j] = rw_state_flux[k,j] * (npts - 1)
                        flux_evol[iblock]['iter_start'][k,j] = start
                        flux_evol[iblock]['iter_stop'][k,j] = stop

                color_prob_evol[iblock] = rw_color_probs
                state_prob_evol[iblock] = rw_state_probs[:-1]
                bin_prob_evol[iblock] = rw_bin_probs

                # Save the transition matrix if desired.
                if self.save_transition_matrices:
                    iter_group = self.output_file['iterations']\
                            .create_group('iter_{:08d}'.format(stop-1))

                    # transition_matrix should already have the 'NaN' elements
                    # set to zero.
                    # Convert to SciPy's sparse coordinate matrix format 
                    sparse_transmat = sp.coo_matrix(transition_matrix)
                     
                    row = sparse_transmat.row
                    col = sparse_transmat.col
                    k   = sparse_transmat.data

                    iter_group.create_dataset('rows', data=row, compression=9)
                    iter_group.create_dataset('cols', data=col, compression=9)
                    iter_group.create_dataset('k',    data=k,   compression=9)
                    iter_group.attrs.create('iter_start', start)
                    iter_group.attrs.create('iter_stop', stop)
                                            
            # Save the data sets
            ds_flux_evol = self.output_file.create_dataset('conditional_flux_evolution', data=flux_evol, shuffle=True, compression=9)
            ds_state_prob_evol = self.output_file.create_dataset('state_prob_evolution', data=state_prob_evol, compression=9)
            ds_color_prob_evol = self.output_file.create_dataset('color_prob_evolution', data=color_prob_evol, compression=9)
            ds_bin_prob_evol = self.output_file.create_dataset('bin_prob_evolution', data=bin_prob_evol, compression=9)
            ds_state_labels = self.output_file.create_dataset('state_labels', data=state_labels)


if __name__ == '__main__':
    WMultiReweightTool().main()
