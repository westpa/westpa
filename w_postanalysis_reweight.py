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
import networkx as nx
import scipy.sparse as sp
import h5py
from h5py import h5s
import sys
import traceback

import westpa
from west.data_manager import weight_dtype, n_iter_dtype
from westtools import (WESTTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)
from westpa import h5io


def normalize(m):
    nm = m.copy()

    row_sum = m.sum(1)
    ii = np.nonzero(row_sum)[0]
    nm[ii,:] = m[ii,:] / row_sum[ii][:, np.newaxis]

    return nm


def steadystate_solve(K):
    # Reformulate K to remove sink/source states
    G = nx.from_numpy_matrix(K, create_using=nx.DiGraph())
    components = nx.strongly_connected_components(G)[0]

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


def accumulate_statistics(h5file, start_iter, stop_iter, nbins, total_fluxes=None, total_obs=None):
    if total_fluxes is None:
        assert total_obs is None
        total_fluxes = np.zeros((nbins, nbins), weight_dtype)
        total_obs = np.zeros((nbins, nbins), np.int64)

    for iiter in xrange(start_iter, stop_iter):
        iter_grp = h5file['iterations']['iter_{:08d}'.format(iiter)]

        rows = iter_grp['rows'][...]
        cols = iter_grp['cols'][...]
        obs = iter_grp['obs'][...]
        flux = iter_grp['flux'][...]

        total_fluxes += sp.coo_matrix((flux, (rows, cols)), shape=(nbins, nbins)).todense()
        total_obs += sp.coo_matrix((obs, (rows, cols)), shape=(nbins, nbins)).todense()

    total_pop = np.sum(h5file['bin_populations'][start_iter:stop_iter, :], axis=0)

    return total_fluxes, total_obs, total_pop


def reweight(h5file, start, stop, nstates, nbins, state_labels, state_map, nfbins, total_fluxes=None, total_obs=None):

    total_fluxes, total_obs, total_pop = accumulate_statistics(h5file, start, stop, nfbins, total_fluxes, total_obs)

    flux_matrix = total_fluxes.copy()
    flux_matrix[total_obs < self.obs_threshold] = 0.0
    transition_matrix = normalize(flux_matrix)

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


def _calc_state_flux(trans_matrix, index1, index2, bin_probs, bin_last_state_map, bin_state_map, nstates):
    state_flux = np.zeros((nstates, nstates), np.float64)
    
    n_trans = index1.shape[0]
    for k in xrange(n_trans):
        ii = bin_last_state_map[index1[k]]
        jj = bin_state_map[index2[k]]

        if jj != nstates:
            state_flux[ii, jj] += trans_matrix[k] * bin_probs[index1[k]]

    return state_flux

try:
    import numba as nb
    calc_state_flux = nb.jit(nb.f8[:,:](nb.f8[:], nb.i8[:], nb.i8[:], nb.f8[:], nb.u2[:], nb.u2[:], nb.i8), nopython=False)(_calc_state_flux)
except:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
    print(''.join(line for line in lines))
    print('Warning: Numba module not present; using slower pure python code to calculate state fluxes')
    calc_state_flux = _calc_state_flux


log = logging.getLogger('w_postanalysis_reweighting')

from westtools.dtypes import iter_block_ci_dtype as ci_dtype


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
        iogroup.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')

        iogroup.add_argument('-k', '--kinetics', default='flux_matrices.h5',
                            help='''Per-iteration flux matrices calculated by w_postanalysis_matrix 
                            (default: %(default)s).''')
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

        cogroup.add_argument('--obs-threshold', type=int, default=1,
                             help='''The minimum number of observed transitions between two states i and j necessary to include
                             fluxes in the reweighting estimate''')
        
    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_filename, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.assignments_file = h5io.WESTPAH5File(self.assignments_filename, 'r')#, driver='core', backing_store=False)
        self.kinetics_file = h5io.WESTPAH5File(self.kinetics_filename, 'r')#, driver='core', backing_store=False)
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments data do not span the requested iterations')

        if not self.iter_range.check_data_iter_range_least(self.kinetics_file):
            raise ValueError('kinetics data do not span the requested iterations')

    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args, default_iter_step=None)
        if self.iter_range.iter_step is None:
            #use about 10 blocks by default
            self.iter_range.iter_step = max(1, (self.iter_range.iter_stop - self.iter_range.iter_start) // 10)
        
        self.output_filename = args.output
        self.assignments_filename = args.assignments
        self.kinetics_filename = args.kinetics
                
        self.evolution_mode = args.evolution_mode
        self.evol_window_frac = args.window_frac
        if self.evol_window_frac <= 0 or self.evol_window_frac > 1:
            raise ValueError('Parameter error -- fractional window defined by --window-frac must be in (0,1]')
        self.obs_threshold = args.obs_threshold



    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            nstates = self.assignments_file.attrs['nstates']
            nbins = self.assignments_file.attrs['nbins']
            state_labels = self.assignments_file['state_labels'][...]
            state_map = self.assignments_file['state_map'][...]
            nfbins = self.kinetics_file.attrs['nrows']

            assert nstates == len(state_labels)
            assert nfbins == nbins * nstates

            start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

            start_pts = range(start_iter, stop_iter, step_iter)
            flux_evol = np.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
            color_prob_evol = np.zeros((len(start_pts), nstates))
            state_prob_evol = np.zeros((len(start_pts), nstates))
            bin_prob_evol = np.zeros((len(start_pts), nfbins))
            pi.new_operation('Calculating flux evolution', len(start_pts))

            if self.evolution_mode == 'cumulative' and self.evol_window_frac == 1.0:
                print('Using fast streaming accumulation')

                total_fluxes = np.zeros((nfbins, nfbins), weight_dtype)
                total_obs = np.zeros((nfbins, nfbins), np.int64)

                for iblock, start in enumerate(start_pts):
                    pi.progress += 1
                    stop = min(start + step_iter, stop_iter)

                    params = dict(start=start, stop=stop, nstates=nstates, nbins=nbins,
                                  state_labels=state_labels, state_map=state_map, nfbins=nfbins,
                                  total_fluxes=total_fluxes, total_obs=total_obs,
                                  h5file=self.kinetics_file)

                    rw_state_flux, rw_color_probs, rw_state_probs, rw_bin_probs, rw_bin_flux = reweight(**params)
                    for k in xrange(nstates):
                        for j in xrange(nstates):
                            flux_evol[iblock]['expected'][k,j] = rw_state_flux[k,j]
                            flux_evol[iblock]['iter_start'][k,j] = start
                            flux_evol[iblock]['iter_stop'][k,j] = stop

                    color_prob_evol[iblock] = rw_color_probs
                    state_prob_evol[iblock] = rw_state_probs[:-1]
                    bin_prob_evol[iblock] = rw_bin_probs


            else:
                for iblock, start in enumerate(start_pts):
                    pi.progress += 1
                    stop = min(start + step_iter, stop_iter)
                    if self.evolution_mode == 'cumulative':
                        windowsize = int(self.evol_window_frac * (stop - start_iter))
                        block_start = max(start_iter, stop - windowsize)
                    else:   # self.evolution_mode == 'blocked'
                        block_start = start

                    params = dict(start=block_start, stop=stop, nstates=nstates, nbins=nbins,
                                  state_labels=state_labels, state_map=state_map, nfbins=nfbins,
                                  total_fluxes=None, total_obs=None,
                                  h5file=self.kinetics_file)

                    rw_state_flux, rw_color_probs, rw_state_probs, rw_bin_probs, rw_bin_flux = reweight(**params)
                    for k in xrange(nstates):
                        for j in xrange(nstates):
                            flux_evol[iblock]['expected'][k,j] = rw_state_flux[k,j]
                            flux_evol[iblock]['iter_start'][k,j] = start
                            flux_evol[iblock]['iter_stop'][k,j] = stop

                    color_prob_evol[iblock] = rw_color_probs
                    state_prob_evol[iblock] = rw_state_probs[:-1]
                    bin_prob_evol[iblock] = rw_bin_probs


            ds_flux_evol = self.output_file.create_dataset('conditional_flux_evolution', data=flux_evol, shuffle=True, compression=9)
            ds_state_prob_evol = self.output_file.create_dataset('state_prob_evolution', data=state_prob_evol, compression=9)
            ds_color_prob_evol = self.output_file.create_dataset('color_prob_evolution', data=color_prob_evol, compression=9)
            ds_bin_prob_evol = self.output_file.create_dataset('bin_prob_evolution', data=bin_prob_evol, compression=9)
            ds_state_labels = self.output_file.create_dataset('state_labels', data=state_labels)


if __name__ == '__main__':
    WPostAnalysisReweightTool().main()
