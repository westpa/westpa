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
import numpy
import scipy.sparse as sp
from scipy.sparse import csgraph
import h5py

from collections import Counter
from postanalysis import reweight_for_c

import westpa
from west.data_manager import weight_dtype, n_iter_dtype
#from westtools import (WESTTool, WESTParallelTool, WESTDataReader, IterRangeSelection,
from westtools import (WESTMasterCommand, WESTParallelTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)

from westpa.kintool import WESTKinAvg, AverageCommands
from westpa import h5io
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

log = logging.getLogger('westtools.w_postanalysis_reweight')

import mclib

from mclib import mcbs_correltime, mcbs_ci_correl_rw, _1D_simple_eval_block, _2D_simple_eval_block

def _2D_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, **extra):
    # Our rate estimator is a little more complex, so we've defined a custom evaluation block for it,
    # instead of just using the block evalutors that we've imported.
    results = []
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            if istate == jstate: continue
            kwargs = { 'istate' : istate, 'jstate': jstate }
            # Ergo, we need to send in... nbins, state_map, return_flux, return_states, return_color.  By default, we always return rates
            kwargs = dict(istate=istate, jstate=jstate, nstates=nstates, nbins=extra['nbins'], state_map=extra['state_map'], return_flux=extra['return_flux'], return_states=extra['return_states'], return_color=extra['return_color'])

            dataset = {'dataset': data_input['dataset'] }
            ci_res = mcbs_ci_correl_rw(dataset,estimator=reweight_for_c,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=(lambda x: x[0]), pre_calculated=data_input['pre_calculated'][:,istate,jstate][numpy.isfinite(data_input['pre_calculated'][:,istate,jstate])], do_correl=do_correl, **kwargs)
            results.append((name, iblock, istate, jstate, (start,stop) + ci_res))

    return results

def _eval_block(iblock, start, stop, nstates, nbins, mcbs_alpha, mcbs_nsets, mcbs_acalpha, rates, flux_c, pop_c, state_map, correl, **kwargs):
    results = [[],[],[],[]]
    # results are conditional fluxes, rates
    for istate in xrange(nstates):
        dataset = { 'indices' : np.array(range(start-1, stop-1), dtype=np.uint16) }
        # We need a jstate for the function, but it isn't really important.
        kwargs.update({ 'istate' : istate, 'jstate': 1, 'nstates' : nstates, 'nbins' : nbins , 'state_map': state_map, 'return_flux': False, 'return_states': True})
        ci_res = mcbs_ci_correl_rw(dataset, estimator=reweight_for_c,
                                alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                subsample=(lambda x: x[0]), pre_calculated=pop_c[:,istate], correl=correl, **kwargs)
        results[0].append((iblock, istate, (start,stop) + ci_res))

        dataset = { 'indices' : np.array(range(start-1, stop-1), dtype=np.uint16) }
        # We need a jstate for the function, but it isn't really important.
        kwargs.update({ 'istate' : istate, 'jstate': 1, 'nstates' : nstates, 'nbins' : nbins , 'state_map': state_map, 'return_flux': False, 'return_states': False, 'return_color': True})
        ci_res = mcbs_ci_correl_rw(dataset, estimator=reweight_for_c,
                                alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                subsample=(lambda x: x[0]), pre_calculated=pop_c[:,istate], correl=correl, **kwargs)
        results[1].append((iblock, istate, (start,stop) + ci_res))

        for jstate in xrange(nstates):
            # Fluxes!
            if istate == jstate: continue
            dataset = { 'indices' : np.array(range(start-1, stop-1), dtype=np.uint16) }
            kwargs.update({ 'istate' : istate, 'jstate': jstate, 'nstates' : nstates, 'nbins' : nbins , 'state_map': state_map, 'return_flux': True, 'return_states': False})
            ci_res = mcbs_ci_correl_rw(dataset, estimator=reweight_for_c,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=(lambda x: x[0]), pre_calculated=flux_c[:,istate,jstate], correl=correl, **kwargs)
            results[2].append((iblock, istate, jstate, (start,stop) + ci_res))

            # Rates!
            if istate == jstate: continue
            dataset = { 'indices' : np.array(range(start-1, stop-1), dtype=np.uint16) }
            kwargs.update({ 'istate' : istate, 'jstate': jstate, 'nstates' : nstates, 'nbins' : nbins , 'state_map': state_map, 'return_flux': False, 'return_states': False})
            ci_res = mcbs_ci_correl_rw(dataset, estimator=reweight_for_c,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=(lambda x: x[0]), pre_calculated=rates[:,istate,jstate], correl=correl, **kwargs)
            results[3].append((iblock, istate, jstate, (start,stop) + ci_res))

    return results

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


def accumulate_statistics(h5file, start_iter, stop_iter, nbins, total_fluxes=None, total_obs=None):
    if total_fluxes is None:
        assert total_obs is None
        total_fluxes = np.zeros((nbins, nbins), weight_dtype)
        total_obs = np.zeros((nbins, nbins), np.int64)

    rows = []
    cols = []
    obs = []
    flux = []

    for iiter in xrange(start_iter, stop_iter):
        iter_grp = h5file['iterations']['iter_{:08d}'.format(iiter)]

        rows.append(iter_grp['rows'][...])
        cols.append(iter_grp['cols'][...])
        obs.append(iter_grp['obs'][...])
        flux.append(iter_grp['flux'][...])

    rows, cols, obs, flux = map(np.hstack, [rows, cols, obs, flux])

    total_fluxes += sp.coo_matrix((flux, (rows, cols)), shape=(nbins, nbins)).todense()
    total_obs += sp.coo_matrix((obs, (rows, cols)), shape=(nbins, nbins)).todense()

    total_pop = np.sum(h5file['bin_populations'][start_iter:stop_iter, :], axis=0)

    return total_fluxes, total_obs, total_pop

def accumulate_statistics_list(h5file, iterations, nbins, total_fluxes=None, total_obs=None):
    if total_fluxes is None:
        assert total_obs is None
        total_fluxes = np.zeros((nbins, nbins), weight_dtype)
        total_obs = np.zeros((nbins, nbins), np.int64)

    rows = []
    cols = []
    obs = []
    flux = []

    for iiter in xrange(iterations):
        iter_grp = h5file['iterations']['iter_{:08d}'.format(iiter)]

        rows.append(iter_grp['rows'][...])
        cols.append(iter_grp['cols'][...])
        obs.append(iter_grp['obs'][...])
        flux.append(iter_grp['flux'][...])

    rows, cols, obs, flux = map(np.hstack, [rows, cols, obs, flux])

    total_fluxes += sp.coo_matrix((flux, (rows, cols)), shape=(nbins, nbins)).todense()
    total_obs += sp.coo_matrix((obs, (rows, cols)), shape=(nbins, nbins)).todense()

    total_pop = np.sum(h5file['bin_populations'][start_iter:stop_iter, :], axis=0)

    return total_fluxes, total_obs, total_pop


def reweight(h5file, start, stop, nstates, nbins, state_labels, state_map, nfbins, obs_threshold=1, total_fluxes=None, total_obs=None):

    # Instead of pulling in start and stop, we'll pull in a list of indices.
    # This way, it should support the bootstrap.
    total_fluxes, total_obs, total_pop = accumulate_statistics(h5file, start, stop, nfbins, total_fluxes, total_obs)

    flux_matrix = total_fluxes.copy()
    flux_matrix[total_obs < obs_threshold] = 0.0
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


def calc_state_flux(trans_matrix, index1, index2, bin_probs, bin_last_state_map, bin_state_map, nstates):
    state_flux = np.zeros((nstates, nstates), np.float64)
    
    n_trans = index1.shape[0]
    for k in xrange(n_trans):
        ii = bin_last_state_map[index1[k]]
        jj = bin_state_map[index2[k]]

        if jj != nstates:
            state_flux[ii, jj] += trans_matrix[k] * bin_probs[index1[k]]

    return state_flux




class WPostAnalysisReweightTool(WESTParallelTool):
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

        self.mcbs_enable = None
        self.mcbs_alpha = None
        self.mcbs_acalpha = None
        self.mcbs_nsets = None
        
    def stamp_mcbs_info(self, dataset):
        dataset.attrs['mcbs_alpha'] = self.mcbs_alpha
        dataset.attrs['mcbs_acalpha'] = self.mcbs_acalpha
        dataset.attrs['mcbs_nsets'] = self.mcbs_nsets
        
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

        cgroup = parser.add_argument_group('confidence interval calculation options')
        cgroup.add_argument('--bootstrap', dest='bootstrap', action='store_const', const=True,
                             help='''Enable the use of Monte Carlo Block Bootstrapping.''')
        cgroup.add_argument('--disable-correl', '-dc', dest='correl', action='store_const', const=False,
                             help='''Disable the correlation analysis.''')
        cgroup.add_argument('--alpha', type=float, default=0.05, 
                             help='''Calculate a (1-ALPHA) confidence interval'
                             (default: %(default)s)''')
        cgroup.add_argument('--autocorrel-alpha', type=float, dest='acalpha', metavar='ACALPHA',
                             help='''Evaluate autocorrelation to (1-ACALPHA) significance.
                             Note that too small an ACALPHA will result in failure to detect autocorrelation
                             in a noisy flux signal. (Default: same as ALPHA.)''')
        cgroup.add_argument('--nsets', type=int, default=1000,
                             help='''Use NSETS samples for bootstrapping (default: chosen based on ALPHA)''')
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

        agroup = parser.add_argument_group('other options')
        agroup.add_argument('--config-from-file', dest='config_from_file', action='store_true', 
                            help='''Load bins/macrostates from a scheme specified in west.cfg.''')
        agroup.add_argument('--scheme-name', dest='scheme',
                            help='''Name of scheme specified in west.cfg.''')
        
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

        self.mcbs_enable = args.bootstrap
        self.mcbs_alpha = args.alpha
        self.mcbs_acalpha = args.acalpha if args.acalpha else self.mcbs_alpha
        self.mcbs_nsets = args.nsets if args.nsets else mclib.get_bssize(self.mcbs_alpha)
        self.correl = args.correl
                
        self.evolution_mode = args.evolution_mode
        self.evol_window_frac = args.window_frac
        if self.evol_window_frac <= 0 or self.evol_window_frac > 1:
            raise ValueError('Parameter error -- fractional window defined by --window-frac must be in (0,1]')
        self.obs_threshold = args.obs_threshold

        if args.config_from_file:
            if not args.scheme:
                raise ValueError('A scheme must be specified.')
            else:
                self.load_config_from_west(args.scheme)

    def load_config_from_west(self, scheme):
        try:
            config = westpa.rc.config['west']['w_ipython']
        except:
            raise ValueError('There is no configuration file specified.')
        import os
        path = os.path.join(os.getcwd(), config['directory'], scheme)
        try:
            os.mkdir(config['directory'])
            os.mkdir(path)
        except:
            pass
        self.output_filename = os.path.join(path, 'kinrw.h5')
        self.kinetics_filename = os.path.join(path, 'flux_matrices.h5')
        self.assignments_filename = os.path.join(path, 'assign.h5')
        w_kinavg_config = { 'mcbs_alpha': 0.05, 'mcbs_nsets': 1000, 'evolution': 'cumulative', 'evol_window_frac': 1, 'step_iter': 1, 'bootstrap': True , 'correl': True, 'obs_threshold': 1 }
        try:
            w_kinavg_config.update(config['w_kinavg'])
        except:
            pass
        try:
            w_kinavg_config.update(config['analysis_schemes'][scheme]['w_kinavg'])
        except:
            pass
        self.mcbs_alpha = w_kinavg_config['mcbs_alpha']
        # Probably problematic, as we should allow this option itself, but there it is for now.
        self.mcbs_acalpha = self.mcbs_alpha
        self.mcbs_nsets = w_kinavg_config['mcbs_nsets']
        self.evolution_mode = w_kinavg_config['evolution']
        self.evol_window_frac = w_kinavg_config['evol_window_frac']
        self.iter_range.iter_step = w_kinavg_config['step_iter']
        self.mcbs_enable = w_kinavg_config['bootstrap']
        self.correl = w_kinavg_config['correl']
        self.obs_threshold = w_kinavg_config['obs_threshold']




    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
                #flux_evol[:]['expected'][k,j] *= (npts - 1)


            if True:
                for iblock, start in enumerate(start_pts):
                    pi.progress += 1
                    
                    stop = min(start + step_iter, stop_iter)
                    if self.evolution_mode == 'cumulative':
                        windowsize = max(1, int(self.evol_window_frac * (stop - start_iter)))
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
                            # Normalize such that we report the flux per tau (tau being the weighted ensemble iteration)
                            # npts DOES NOT always includes a 0th time point; must be a way to gracefully handle this?
                            flux_evol[iblock]['expected'][k,j] = rw_state_flux[k,j] * (npts - 1)
                            flux_evol[iblock]['iter_start'][k,j] = start
                            flux_evol[iblock]['iter_stop'][k,j] = stop
                            rate_evol[iblock]['expected'][k,j] = (rw_state_flux[k,j] * (npts - 1)) / rw_color_probs[k]
                            if rw_color_probs[k] == 0.0  or flux_evol[iblock]['expected'][k,j] == 0.0:
                                rate_evol[iblock]['expected'][k,j] = 0
                            rate_evol[iblock]['iter_start'][k,j] = start
                            rate_evol[iblock]['iter_stop'][k,j] = stop

                    color_prob_evol[iblock] = rw_color_probs
                    state_prob_evol[iblock] = rw_state_probs[:-1]
                    bin_prob_evol[iblock] = rw_bin_probs

            if self.mcbs_enable == True:
                futures = []
                all_items = np.arange(1,len(start_pts)+1)
                bootstrap_length = 0.5*(len(start_pts)*(len(start_pts)+1)) - np.delete(all_items, np.arange(1, len(start_pts)+1, step_iter))
                pi.new_operation('Calculating Monte Carlo Bootstrap', bootstrap_length[0])
                rows = []
                cols = []
                obs = []
                flux = []
                insert = []

                for iiter in xrange(start_iter, stop_iter):
                    iter_grp = self.kinetics_file['iterations']['iter_{:08d}'.format(iiter)]

                    rows.append(iter_grp['rows'][...])
                    cols.append(iter_grp['cols'][...])
                    obs.append(iter_grp['obs'][...])
                    flux.append(iter_grp['flux'][...])
                    if iiter != start_iter:
                        insert.append(iter_grp['rows'][...].shape[0] + insert[-1])
                    else:
                        insert.append(iter_grp['rows'][...].shape[0])
                rows = np.concatenate(rows)
                cols = np.concatenate(cols)
                obs = np.concatenate(obs)
                assert insert[-1] == len(rows)
                flux = np.concatenate(flux)
                ins = []
                ins.append(0)
                ins += insert
                insert = np.array(ins, dtype=np.intc)

        
                for iblock, start in enumerate(start_pts):
                    stop = min(start+step_iter, stop_iter)
                    if self.evolution_mode == 'cumulative':
                        windowsize = max(1, int(self.evol_window_frac * (stop - start_iter)))
                        block_start = max(start_iter, stop - windowsize)
                    else:   # self.evolution_mode == 'blocked'
                        block_start = start
                    # Why flux_c?
                    # This is one where we've already calculated the fluxes.  I can probably set it to be different, but.
                    future = self.work_manager.submit(_eval_block, kwargs=dict(iblock=iblock, start=block_start, stop=stop,
                                                                               nstates=nstates, nbins=nbins, nfbins=nfbins,
                                                                               rows=rows, cols=cols, obs=obs, flux=flux, insert=insert,
                                                                               state_labels=state_labels, 
                                                                               bin_last_state_map=np.tile(np.arange(nstates, dtype=np.int), nbins), 
                                                                               bin_state_map=np.repeat(state_map[:-1], nstates),
                                                                               state_map=state_map,
                                                                               rates=rate_evol[block_start-1:stop,:,:]['expected'],
                                                                               flux_c=flux_evol[block_start-1:stop,:,:]['expected'],
                                                                               pop_c=state_prob_evol[block_start-1:stop,:],
                                                                               mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                                                               mcbs_acalpha=self.mcbs_acalpha,
                                                                               correl=self.correl))
                    futures.append(future)
                
                for future in self.work_manager.as_completed(futures):
                    state_results, color_results, condflux_results, rate_results = future.get_result(discard=True)

                    for result in state_results:
                        iblock,istate,ci_result = result
                        state_prob_bootstrap[iblock, istate] = ci_result

                    for result in color_results:
                        iblock,istate,ci_result = result
                        color_prob_bootstrap[iblock, istate] = ci_result
                        
                    for result in condflux_results:
                        iblock,istate,jstate,ci_result = result
                        flux_evol_bootstrap[iblock, istate, jstate] = ci_result
                        # Normalisation for the pcoord length, to ensure the results are per tau, not per subtau.
                        # There should be a way to gracefully handle simulations without a 0th timepoint.
                        flux_evol_bootstrap[iblock, istate, jstate]['expected'] = ci_result[2] * (npts - 1)
                        flux_evol_bootstrap[iblock, istate, jstate]['ci_lbound'] = ci_result[3] * (npts - 1)
                        flux_evol_bootstrap[iblock, istate, jstate]['ci_ubound'] = ci_result[4] * (npts - 1)
                    
                    for result in rate_results:
                        iblock, istate, jstate, ci_result = result 
                        rate_evol_bootstrap[iblock, istate, jstate] = ci_result
                        rate_evol_bootstrap[iblock, istate, jstate]['expected'] = ci_result[2] * (npts - 1)
                        rate_evol_bootstrap[iblock, istate, jstate]['ci_lbound'] = ci_result[3] * (npts - 1)
                        rate_evol_bootstrap[iblock, istate, jstate]['ci_ubound'] = ci_result[4] * (npts - 1)


                    pi.progress += iblock / step_iter
            else:
                flux_evol_bootstrap = flux_evol
                rate_evol_bootstrap = rate_evol
                state_prob_bootstrap = state_prob_evol
                color_prob_bootstrap = color_prob_evol



            ds_flux_evol = self.output_file.create_dataset('conditional_flux_evolution', data=flux_evol_bootstrap, shuffle=True, compression=9)
            ds_flux_evol = self.output_file.create_dataset('rate_evolution', data=rate_evol_bootstrap, shuffle=True, compression=9)
            ds_state_prob_evol = self.output_file.create_dataset('state_pop_evolution', data=state_prob_bootstrap, compression=9)
            ds_color_prob_evol = self.output_file.create_dataset('color_prob_evolution', data=color_prob_bootstrap, compression=9)
            ds_bin_prob_evol = self.output_file.create_dataset('bin_prob_evolution', data=bin_prob_evol, compression=9)
            ds_state_labels = self.output_file.create_dataset('state_labels', data=state_labels)

            if self.mcbs_enable == True:
                for ds in (ds_flux_evol, ds_state_prob_evol, ds_color_prob_evol, ds_bin_prob_evol):
                    self.stamp_mcbs_info(ds)

class RWReweight(AverageCommands):
    subcommand = 'reweight'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'direct.h5'
    description = '''\ Blah!'''

    def generate_reweight_data(self):

        # We're initializing the various datasets...
        if True:
            self.open_files()
            self.open_assignments()

            # This is the reweighting specific code.  We want to generate the rates like we normally would...
            # Eventually, we want all our datasets to be like this.

            self.nfbins = self.kinetics_file.attrs['nrows']
            self.npts = self.kinetics_file.attrs['npts']

            assert self.nstates == len(self.state_labels)
            assert self.nfbins == self.nbins * self.nstates

            start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

            start_pts = range(start_iter, stop_iter, step_iter)

            rates = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates, self.nstates), dtype=ci_dtype))
            flux = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates, self.nstates), dtype=ci_dtype))
            color_prob = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates), dtype=ci_dtype))
            state_prob = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates), dtype=ci_dtype))
            bin_prob = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nfbins), dtype=ci_dtype))

            # We're pre-generating the datasets for use with the estimators.
            # This means running the reweighting code through all the iterations selected, as our mclib portion requires
            # that we have this dataset pre-generated.

            if True:

                total_fluxes = np.zeros((self.nfbins, self.nfbins), weight_dtype)
                total_obs = np.zeros((self.nfbins, self.nfbins), np.int64)

                for iblock, start in enumerate(start_pts):
                    pi.progress += 1
                    stop = min(start + step_iter, stop_iter)

                    params = dict(start=start, stop=stop, nstates=self.nstates, nbins=self.nbins,
                                  state_labels=self.state_labels, state_map=self.state_map, nfbins=self.nfbins,
                                  total_fluxes=total_fluxes, total_obs=total_obs,
                                  h5file=self.kinetics_file, obs_threshold=self.obs_threshold)

                    rw_state_flux, rw_color_probs, rw_state_probs, rw_bin_probs, rw_bin_flux = reweight(**params)

                    # Now, set the data.
                    for k in xrange(nstates):
                        for j in xrange(nstates):
                            # We need to take the iteration into account, as well...
                            # Normalize such that we report the flux per tau (tau being the weighted ensemble iteration)
                            # npts MAY include a 0th time point

                            rates.data[iblock,k,j] = (rw_state_flux[k,j] * (npts - 1)) / rw_color_probs[k]
                            if rw_color_probs[k] == 0.0  or rw_state_flux[k,j] == 0.0:
                                rates.data = 0
                            flux.data[iblock,k,j] = rw_state_flux[k,j] * (npts - 1)
                            color_prob.data[iblock] = rw_color_probs
                            state_prob.data[iblock] = rw_state_probs[:-1]
                            bin_prob.data[iblock] = rw_bin_probs

                # Now we're setting the data, here...

                self.rates = rates
                self.flux = flux
                self.color_prob = color_prob
                self.state_prob = state_prob
                self.bin_prob = bin_prob

    def postanalysis_reweight(self):
        pi = self.progress.indicator

        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

        start_pts = range(start_iter, stop_iter, step_iter)


        indices = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts)), dtype=ci_dtype))
        indices.data = np.array(range(start_iter-1, stop_iter-1), dtype=np.uint16)

        submit_kwargs = dict(pi=pi, nstates=self.nstates, start_iter=self.start_iter, stop_iter=self.stop_iter, 
                             step_iter=self.step_iter, nbins=self.nbins, state_map=self.state_map)

        # Due to the way our estimator is written, we're using the same dataset every time.  We're just returning different values.
        submit_kwargs['dataset'] = {'dataset': indices}


        # Calculate averages for the simulation, then report, if necessary.
        with pi:

            submit_kwargs.update(dict(return_flux=False, return_states=False, return_color=False, pre_calculated=self.rates))
            avg_rates = self.run_calculation(eval_block=_2D_eval_block, name='Rate Evolution', dim=2, do_averages=True, **submit_kwargs)
            self.output_file.replace_dataset('avg_rates', data=avg_rates[1])

            submit_kwargs.update(dict(return_flux=True, return_states=False, return_color=False, pre_calculated=self.flux))
            avg_conditional_fluxes = self.run_calculation(eval_block=_2D_eval_block, name='Conditional Flux Evolution', dim=2, do_averages=True, **submit_kwargs)
            self.output_file.replace_dataset('avg_conditional_fluxes', data=avg_conditional_fluxes[1])

        # Now, print them!
        pi.clear()

        # We've returned an average, but it still exists in a timeslice format.  So we need to return the 'last' value.
        self.print_averages(avg_conditional_fluxes[1], '\nfluxes from state to state:', dim=2)
        self.print_averages(avg_rates[1], '\nrates from state to state:', dim=2)

        # Do a bootstrap evolution.
        pi.clear()
        with pi:
            submit_kwargs.update(dict(return_flux=False, return_states=False, return_color=False, pre_calculated=self.rates))
            rate_evol = self.run_calculation(eval_block=_2D_eval_block, name='Rate Evolution', dim=2, **submit_kwargs)
            self.output_file.replace_dataset('rate_evolution', data=rate_evol, shuffle=True, compression=9)
            pi.clear()

            submit_kwargs.update(dict(return_flux=True, return_states=False, return_color=False, pre_calculated=self.flux))
            rate_evol = self.run_calculation(eval_block=_2D_eval_block, name='Conditional Flux Evolution', dim=2, **submit_kwargs)
            self.output_file.replace_dataset('conditional_flux_evolution', data=rate_evol, shuffle=True, compression=9)
            pi.clear()

    def go(self):
        self.w_postanalysis_reweight()


class RWStateProbs(AverageCommands):
    subcommand = 'reweight'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'direct.h5'
    description = '''\ Blah!'''
    def w_postanalysis_stateprobs():
        pass

class RWFlux(AverageCommands):
    subcommand = 'flux'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'direct.h5'
    description = '''\ Blah!'''
    def w_postanalysis_stateprobs():
        pass

class RWAll(RWStateProbs, RWReweight, RWFlux):
    subcommand = 'all'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'direct.h5'

    def go(self):
        # One minor issue; as this stands now, since it's inheriting from all the other classes, it needs
        # a kinetics file to instantiate the other attributes.  We'll need to modify how the loading works, there.
        self.w_postanalysis_matrix()
        self.w_postanalysis_reweight()
        self.w_postanalysis_stateprobs()

# Just a convenience class to average the observables.
class RWAverage(RWStateProbs, RWReweight):
    subcommand = 'average'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'direct.h5'

    def go(self):
        self.w_postanalysis_reweight()
        self.w_postanalysis_stateprobs()

class WReweight(WESTMasterCommand, WESTParallelTool):
    prog='w_reweight'
    #subcommands = [AvgTraceSubcommand,AvgMatrixSubcommand]
    subcommands = [RWFlux, RWReweight, RWStateProbs, RWAll, RWAverage]
    subparsers_title = 'direct kinetics analysis schemes'

if __name__ == '__main__':
    WReweight().main()
