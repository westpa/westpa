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
from westtools import (WESTTool, WESTDataReader, IterRangeSelection,
                       ProgressIndicatorComponent)
from westpa import h5io
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

log = logging.getLogger('westtools.w_postanalysis_reweight')

# TODO: add ability to read YAML files (part of westtools core)
# TODO: add in option to actually use the cumulative mean generator


def normalize(m):
    nm = m.copy()

    row_sum = m.sum(1)
    ii = np.nonzero(row_sum)[0]
    nm[ii,:] = m[ii,:] / row_sum[ii][:, np.newaxis]

    return nm


def steadystate_solve(K):
    # Reformulate K to remove sink/source states
    n_components, component_assignments = csgraph.connected_components(K, connection="strong")
    largest_component = Counter(component_assignments).most_common(1)[0][0]
    components = np.where(component_assignments == largest_component)[0]

    ii = np.ix_(components, components)
    K_mod = K[ii]
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

def get_transmat_and_obsmat_sum(fluxmatH5_list, start_iter, stop_iter, nbins, 
                                recycling_bin_list=None):
    '''
    Find the sum of transition matrices from a set of simulations, for a given 
    set of timepoints. 

    Return: 
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

def get_average_transition_mat(fluxmatH5_list, start_iter, stop_iter, nbins, 
                               recycling_bin_list=None): 
    '''
    Find an average transition matrix from a set of simulations, for a given 
    set of timepoints. Return a numpy matrix representing the average 
    transition matrix, with elements for which no observations were made set
    to ``NaN`` (not a number). 

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
            
def transmat_cumulative_mean_generator(fluxmatH5_list, start_iter, stop_iter, 
                                       step_iter, nbins, 
                                       recycling_bin_list=None):
    '''
    Generator for fast calculations of cumulatively averaged transition matrix 
    estimates.

    Find an average transition matrix from a set of simulations, for a given 
    set of timepoints. Return a numpy matrix representing the average 
    transition matrix, with elements for which no observations were made set
    to ``NaN`` (Not a Number). 

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
    # TODO: remove unnecessary arguments.

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




class WPostAnalysisReweightTool(WESTTool):
    prog ='w_postanalysis_reweight'
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


-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WPostAnalysisReweightTool, self).__init__()
        
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
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
        self.iter_range.include_args['iter_step'] = True
        self.iter_range.add_args(parser)

        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('-y', '--yaml', dest=yamlfilepath, 
                             metavar='YAMLFILE', 
                             default='multi_reweight.yaml',
                             help='''Load options from YAMLFILE. For each 
                             simulation, specify an assignments file and flux
                             matrices file. Search for files in 
                             ['simulations'][SIMNAME]['assignments'] and
                             ['simulations'][SIMNAME]['kinetics']
                             '''
        iogroup.add_argument('-o', '--output', dest='output', default='kinrw.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

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
        cogroup.add_argument('--obs-threshold', type=int, default=1,
                             help='''The minimum number of observed transitions between two states i and j necessary to include
                             fluxes in the reweighting estimate''')
                                                         
    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args, default_iter_step=None)
        if self.iter_range.iter_step is None:
            #use about 10 blocks by default
            self.iter_range.iter_step = max(1, (self.iter_range.iter_stop - self.iter_range.iter_start) // 10)
        
        self.output_filename = args.output
        #self.assignments_filename = args.assignments
        #self.kinetics_filename = args.kinetics
                
        self.evolution_mode = args.evolution_mode
        self.evol_window_frac = args.window_frac
        if self.evol_window_frac <= 0 or self.evol_window_frac > 1:
            raise ValueError('Parameter error -- fractional window defined by --window-frac must be in (0,1]')
        self.obs_threshold = args.obs_threshold

        # Get the list of recycling bins from the YAML file, setting to None
        # if not specified.
        self.recycling_bin_list = [] 
        for simname in self.yamlargdict['simulations'].keys():
            simdict = self.yamlargdict['simulations'][simname]
            if 'recycling_bin_list' in simdict.keys(): 
                self.recycling_bin_list.append(simdict['recycling_bin_list'])
            else:
                self.recycling_bin_list.append(None)
                        

    def open_files(self):
        '''
        Create the output file.  Open input files, including assignments and 
        kinetics files for each simulation specified.  Make the output file 
        available as an attributue of ``self``, ``output_file``.  Similarly,
        make lists of the input files available as attributes, 
        ``assignments_file_list`` and ``kinetics_files_list``.
        '''
        self.output_file = h5io.WESTPAH5File(self.output_filename, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.assignments_file_list = []
        self.kinetics_file_list = []
        self.simname_list = []
        # Open the assignments and kinetics (flux matrices) files for each
        # simulation.  assignments_file_list[i] and kinetcs_file_list[i] 
        # correspond to the same simulation.
        for simname in self.yamldict['simulations'].keys():
            assignments_filename = self.yamldict['simulations'][simname]\
                                                ['assignments']
            assignH5 = h5io.WESTPAH5File(assignments_filename, 'r')
            # Double check that adding this file to the list adds the file
            # itself, and not a pointer to the assignH5 variable.
            self.assignments_file_list.append(assignH5)  
            kinetics_filename = self.yamldict['simulations'][simname]\
                                             ['kinetics']
            kinetH5 = h5io.WESTPAH5File(kinetics_filename, 'r')
            self.kinetics_file_list.append(kinetH5)

            # Keep track of the name associated with each simulation. This tool
            # uses these names later to help give descriptive errors!
            self.simname_list.append(simname)
            # Check that the specified assignment and kinetics (flux matrices)
            # files span the requested iterations.
            if not self.iter_range.check_data_iter_range_least(assignH5):
                raise ValueError('Assignments data {:s} do not span the '
                                 'requested iterations for simulation '
                                 '{:s}'.format(assignments_filename, simname))

            if not self.iter_range.check_data_iter_range_least(kinetH5):
                raise ValueError('Kinetics data {:s} do not span the '
                                 'requested iterations for simulation '
                                 '{:s}'.format(kinetics_filename, simname))

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
            if bincount*statecount != maxrows:
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
             

    def go(self):
        '''Run the main analysis.'''
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            # Get the number of states and the number of bins for the
            # simulations (this should be the same across all the simulations,
            # since this was already checked with 
            # ``check_consistency_of_input_files``).
            nstates = self.assignment_file_list[0].attrs['nstates']
            nbins = self.assignment_file_list[0].attrs['nbins']
            # Get the state labels for each simulation.  These should also be
            # consistent across simulations, but this is up to the user to make
            # sure of.  We could add in a check here later if we decide it
            # makes more sense.
            state_labels = self.assignment_file_list[0].['state_labels'][...]
            # State map *should* also be consistent, but we do not explicitly 
            # enforce that it is.
            state_map = self.assignments_file['state_map'][...]
            # Get nfbins and npts (forced to be consistent by
            # ``check_consistency_of_input_files``).
            nfbins = self.kinetics_file.attrs['nrows']
            npts = self.kinetics_file.attrs['npts']

            assert nstates == len(state_labels)
            assert nfbins == nbins * nstates

            start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

            start_pts = range(start_iter, stop_iter, step_iter)
            flux_evol = np.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
            color_prob_evol = np.zeros((len(start_pts), nstates))
            state_prob_evol = np.zeros((len(start_pts), nstates))
            bin_prob_evol = np.zeros((len(start_pts), nfbins))
            pi.new_operation('Calculating flux evolution', len(start_pts))

            # Set up the generator, if needed; The generator reduces repetitive
            # calculations 
            if self.evolution_mode == 'cumulative' \
                    and self.evol_window_frac == 1.0:
                gen = transmat_cumulative_mean_generator(fluxmatH5_list,
                                                         start_iter, stop_iter,
                                                         step_iter, nbins,
                                                         recycling_bin_list
                                                         ) 
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
                    transition_matrix = get_average_transition_mat(
                                                         fluxmatH5_list
                                                         start_iter, stop_iter,
                                                         nbins,
                                                         recycling_bin_list
                                                                   ) 
                # Apply reweighting procedure 
                rw_state_flux, rw_color_probs, rw_state_probs, \
                        rw_bin_probs, rw_bin_flux = reweight(transition_matrix
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
            # Save the data sets
            ds_flux_evol = self.output_file.create_dataset('conditional_flux_evolution', data=flux_evol, shuffle=True, compression=9)
            ds_state_prob_evol = self.output_file.create_dataset('state_prob_evolution', data=state_prob_evol, compression=9)
            ds_color_prob_evol = self.output_file.create_dataset('color_prob_evolution', data=color_prob_evol, compression=9)
            ds_bin_prob_evol = self.output_file.create_dataset('bin_prob_evolution', data=bin_prob_evol, compression=9)
            ds_state_labels = self.output_file.create_dataset('state_labels', data=state_labels)


if __name__ == '__main__':
    WPostAnalysisReweightTool().main()
