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

# From postanalysis matrix
from westpa.binning import index_dtype
from postanalysis import stats_process

def _2D_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, estimator_kwargs):
    # Our rate estimator is a little more complex, so we've defined a custom evaluation block for it,
    # instead of just using the block evalutors that we've imported.
    results = []
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            if istate == jstate: continue
            kwargs = { 'istate' : istate, 'jstate': jstate }
            # Ergo, we need to send in... nbins, state_map, return_flux, return_states, return_color.  By default, we always return rates
            kwargs = dict(istate=istate, jstate=jstate, nstates=nstates, nbins=estimator_kwargs['nbins'], state_map=estimator_kwargs['state_map'], return_flux=estimator_kwargs['return_flux'], return_states=estimator_kwargs['return_states'], return_color=estimator_kwargs['return_color'])
            kwargs.update(estimator_kwargs)

            #dataset = {'indices': data_input['indices'] }
            dataset = { 'indices' : np.array(range(start-1, stop-1), dtype=np.uint16) }
            pre_calculated = data_input['pre_calculated'][:,istate,jstate]
            #pre_calculated = pre_calculated[np.isfinite(pre_calculated)]
            
            ci_res = mcbs_ci_correl_rw(dataset,estimator=reweight_for_c,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=(lambda x: x[0]), pre_calculated=pre_calculated, do_correl=do_correl, estimator_kwargs=kwargs)
            results.append((name, iblock, istate, jstate, (start,stop) + ci_res))

    return results

def _1D_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, estimator_kwargs):
    # Our rate estimator is a little more complex, so we've defined a custom evaluation block for it,
    # instead of just using the block evalutors that we've imported.
    results = []
    for istate in xrange(nstates):
            kwargs = { 'istate' : istate, 'jstate': istate }
            # Ergo, we need to send in... nbins, state_map, return_flux, return_states, return_color.  By default, we always return rates
            kwargs = dict(istate=istate, jstate=istate, nstates=nstates, nbins=estimator_kwargs['nbins'], state_map=estimator_kwargs['state_map'], return_flux=estimator_kwargs['return_flux'], return_states=estimator_kwargs['return_states'], return_color=estimator_kwargs['return_color'])
            kwargs.update(estimator_kwargs)

            #dataset = {'indices': data_input['indices'] }
            dataset = { 'indices' : np.array(range(start-1, stop-1), dtype=np.uint16) }
            pre_calculated = data_input['pre_calculated'][:,istate]
            #pre_calculated = pre_calculated[np.isfinite(pre_calculated)]
            
            ci_res = mcbs_ci_correl_rw(dataset,estimator=reweight_for_c,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=(lambda x: x[0]), pre_calculated=pre_calculated, do_correl=do_correl, estimator_kwargs=kwargs)
            results.append((name, iblock, istate, (start,stop) + ci_res))

    return results

# From the original postanalysis reweight...

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

# From the original postanalysis matrix

def calc_stats(bin_assignments, weights, fluxes, populations, trans, mask, sampling_frequency):
    fluxes.fill(0.0)
    populations.fill(0.0)
    trans.fill(0)

    stats_process(bin_assignments, weights, fluxes, populations, trans, mask, interval=sampling_frequency)

class RWMatrix(AverageCommands):
    subcommand = 'matrix'
    default_kinetics_file = 'reweight.h5'
    default_output_file = 'reweight.h5'
    help_text = 'create a color-labeled transition matrix from a WESTPA simulation'
    description = '''\
Generate a colored transition matrix from a WE assignment file. The subsequent
analysis requires that the assignments are calculated using only the initial and 
final time points of each trajectory segment. This may require downsampling the
h5file generated by a WE simulation. In the future w_assign may be enhanced to optionally
generate the necessary assignment file from a h5file with intermediate time points.
Additionally, this analysis is currently only valid on simulations performed under
either equilibrium or steady-state conditions without recycling target states.
'''
    
    def __init__(self, parent):
        super(RWMatrix, self).__init__(parent)
        self.progress = ProgressIndicatorComponent()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection() 
        self.output_file = None
        self.assignments_file = None
        #self.default_output_file = 'flux_matrices.h5'
        self.window_size = None
        

    def more_args(self, parser):
        #self.data_reader.add_args(parser)
        #self.iter_range.add_args(parser)
        
        #iogroup = parser.add_argument_group('input/output options')
        #iogroup.add_argument('-a', '--assignments', default='assign.h5',
        #                     help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
        #                        (default: %(default)s).''')

        #iogroup.add_argument('-o', '--output', dest='output', default=self.default_output_file,
        #                     help='''Store results in OUTPUT (default: %(default)s).''')
        cogroup = parser.add_argument_group('calculation options')
        cogroup.add_argument('-s', '--sampling-frequency', 
                             dest='sampling_frequency', 
                             choices=['timepoint','iteration'],
                             default='timepoint',
                             help='''Observe for transition events with a lag
                             time corresponding to either each weighted ensemble
                             iteration (``iteration``) or each sub-iteration 
                             timepoint (``timepoint``) (default: %(default)s).
                             For each pair of bins, store the sum of fluxes
                             in a given iteration (['flux']). Similarly, 
                             for each bin, store the sum of weights in the bin 
                             over all observation points for a given iteration 
                             (['bin_populations']).''' )

                                       
        self.progress.add_args(parser)
        
    def process_more_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)
        if args.config_from_file == False:
            self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
            self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        if args.config_from_file:
            if not args.scheme:
                raise ValueError('A scheme must be specified.')
            else:
                self.load_config_from_west(args.scheme)
        h5io.stamp_creator_data(self.output_file)
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments do not span the requested iterations')
        self.sampling_frequency = args.sampling_frequency

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
        self.output_file = h5io.WESTPAH5File(os.path.join(path, 'flux_matrices.h5'), 'w', creating_program=True)
        self.assignments_file = h5io.WESTPAH5File(os.path.join(path, 'assign.h5'), 'r')
        matrix_config = { 'sampling_frequency': 'timepoint' }
        try:
            matrix_config.update(config['w_postanalysis_matrix'])
        except:
            pass
        try:
            matrix_config.update(config['analysis_schemes'][scheme]['w_postanalysis_matrix'])
        except:
            pass
        self.sampling_frequency = matrix_config['sampling_frequency']


    def w_postanalysis_matrix(self):
        pi = self.progress.indicator
        pi.new_operation('Initializing')
        with pi:

            self.data_reader.open('r')
            nbins = self.assignments_file.attrs['nbins']

            state_labels = self.assignments_file['state_labels'][...]
            state_map = self.assignments_file['state_map'][...]
            nstates = len(state_labels)

            start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
            iter_count = stop_iter - start_iter

            nfbins = nbins * nstates

            flux_shape = (iter_count, nfbins, nfbins)
            pop_shape = (iter_count, nfbins)

            h5io.stamp_iter_range(self.output_file, start_iter, stop_iter)

            bin_populations_ds = self.output_file.create_dataset('bin_populations', shape=pop_shape, dtype=weight_dtype)
            h5io.stamp_iter_range(bin_populations_ds, start_iter, stop_iter)
            h5io.label_axes(bin_populations_ds, ['iteration', 'bin'])

            flux_grp = self.output_file.create_group('iterations')
            self.output_file.attrs['nrows'] = nfbins
            self.output_file.attrs['ncols'] = nfbins

            fluxes = np.empty(flux_shape[1:], weight_dtype)
            populations = np.empty(pop_shape[1:], weight_dtype)
            trans = np.empty(flux_shape[1:], np.int64)

            # Check to make sure this isn't a data set with target states
            #tstates = self.data_reader.data_manager.get_target_states(0)
            #if len(tstates) > 0:
            #    raise ValueError('Postanalysis reweighting analysis does not support WE simulation run under recycling conditions')

            pi.new_operation('Calculating flux matrices', iter_count)
            # Calculate instantaneous statistics
            for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
                # Get data from the main HDF5 file
                iter_group = self.data_reader.get_iter_group(n_iter)
                seg_index = iter_group['seg_index']
                nsegs, npts = iter_group['pcoord'].shape[0:2] 
                weights = seg_index['weight']


                # Get bin and traj. ensemble assignments from the previously-generated assignments file
                assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
                bin_assignments = np.require(self.assignments_file['assignments'][assignment_iiter + np.s_[:nsegs,:npts]],
                                                dtype=index_dtype)

                mask_unknown = np.zeros_like(bin_assignments, dtype=np.uint16)

                macrostate_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
                macrostate_assignments = np.require(self.assignments_file['trajlabels'][macrostate_iiter + np.s_[:nsegs,:npts]],
                                            dtype=index_dtype)

                # Transform bin_assignments to take macrostate membership into account
                bin_assignments  = nstates * bin_assignments + macrostate_assignments

                mask_indx = np.where(macrostate_assignments == nstates)
                mask_unknown[mask_indx] = 1

                # Calculate bin-to-bin fluxes, bin populations and number of obs transitions
                calc_stats(bin_assignments, weights, fluxes, populations, trans, mask_unknown, self.sampling_frequency)

                # Store bin-based kinetics data
                bin_populations_ds[iiter] = populations

                # Setup sparse data structures for flux and obs
                fluxes_sp = sp.coo_matrix(fluxes)
                trans_sp = sp.coo_matrix(trans)

                assert fluxes_sp.nnz == trans_sp.nnz

                flux_iter_grp = flux_grp.create_group('iter_{:08d}'.format(n_iter))
                flux_iter_grp.create_dataset('flux', data=fluxes_sp.data, dtype=weight_dtype)
                flux_iter_grp.create_dataset('obs', data=trans_sp.data, dtype=np.int32)
                flux_iter_grp.create_dataset('rows', data=fluxes_sp.row, dtype=np.int32)
                flux_iter_grp.create_dataset('cols', data=fluxes_sp.col, dtype=np.int32)
                flux_iter_grp.attrs['nrows'] = nfbins
                flux_iter_grp.attrs['ncols'] = nfbins

                # Do a little manual clean-up to prevent memory explosion
                del iter_group, weights, bin_assignments
                del macrostate_assignments

                pi.progress += 1

            # Check and save the number of intermediate time points; this will be used to normalize the
            # flux and kinetics to tau in w_postanalysis_reweight.
            self.output_file.attrs['npts'] = npts

    def go(self):
        self.w_postanalysis_matrix()

class RWReweight(AverageCommands):
    subcommand = 'reweight'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'reweight.h5'
    description = '''\ Blah!'''

    def more_args(self, parser):
        cogroup = parser.add_argument_group('calculation options')
        cogroup.add_argument('--obs-threshold', type=int, default=1,
                             help='''The minimum number of observed transitions between two states i and j necessary to include
                             fluxes in the reweighting estimate''')

    def process_more_args(self, args):
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

    def generate_reweight_data(self):

        # We're initializing the various datasets...
        if True:
            self.open_files()
            self.open_assignments()

            # This is the reweighting specific code.  We want to generate the rates like we normally would...
            # Eventually, we want all our datasets to be like this.

            self.nfbins = self.kinetics_file.attrs['nrows']
            self.npts = self.kinetics_file.attrs['npts']
            self.state_map = self.assignments_file['state_map'][...]
            self.state_labels = self.assignments_file['state_labels'][...]

            assert self.nstates == len(self.state_labels)
            assert self.nfbins == self.nbins * self.nstates

            start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

            start_pts = np.array(range(start_iter, stop_iter, step_iter))

            #rates = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates, self.nstates), dtype=ci_dtype))
            #flux = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates, self.nstates), dtype=ci_dtype))
            #color_prob = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates), dtype=ci_dtype))
            #state_prob = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nstates), dtype=ci_dtype))
            #bin_prob = h5io.IterBlockedDataset.empty_like(np.zeros((len(start_pts), self.nfbins), dtype=ci_dtype))
            rates = np.zeros((len(start_pts), self.nstates, self.nstates), dtype=ci_dtype)
            flux = np.zeros((len(start_pts), self.nstates, self.nstates), dtype=ci_dtype)
            color_prob = np.zeros((len(start_pts), self.nstates), dtype=ci_dtype)
            state_prob = np.zeros((len(start_pts), self.nstates), dtype=ci_dtype)
            bin_prob = np.zeros((len(start_pts), self.nfbins), dtype=ci_dtype)

            # We're pre-generating the datasets for use with the estimators.
            # This means running the reweighting code through all the iterations selected, as our mclib portion requires
            # that we have this dataset pre-generated.

            if True:

                total_fluxes = np.zeros((self.nfbins, self.nfbins), weight_dtype)
                total_obs = np.zeros((self.nfbins, self.nfbins), np.int64)

                for iblock, start in enumerate(start_pts):
                    #pi.progress += 1
                    stop = min(start + step_iter, stop_iter)

                    params = dict(start=start, stop=stop, nstates=self.nstates, nbins=self.nbins,
                                  state_labels=self.state_labels, state_map=self.state_map, nfbins=self.nfbins,
                                  total_fluxes=total_fluxes, total_obs=total_obs,
                                  h5file=self.kinetics_file, obs_threshold=self.obs_threshold)

                    rw_state_flux, rw_color_probs, rw_state_probs, rw_bin_probs, rw_bin_flux = reweight(**params)

                    # Now, set the data.
                    for k in xrange(self.nstates):
                        for j in xrange(self.nstates):
                            # We need to take the iteration into account, as well...
                            # Normalize such that we report the flux per tau (tau being the weighted ensemble iteration)
                            # npts MAY include a 0th time point

                            rates[start:stop,k,j] = (rw_state_flux[k,j] * (self.npts - 1)) / rw_color_probs[k]
                            # We need to handle what the issue is, here!
                            #if rw_color_probs[k] == 0.0  or rw_state_flux[k,j] == 0.0:
                            #    rates.data = 0
                            flux[start:stop,k,j] = rw_state_flux[k,j] * (self.npts - 1)
                            color_prob[start:stop] = rw_color_probs
                            state_prob[start:stop] = rw_state_probs[:-1]
                            bin_prob[start:stop] = rw_bin_probs

                # Now we're setting the data, here...

                self.rates = rates['expected']
                self.flux = flux['expected']
                self.color_prob = color_prob['expected']
                self.state_prob = state_prob['expected']
                self.bin_prob = bin_prob['expected']


                # This is designed to pull in the flux data...
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

                self.rows = rows
                self.cols = cols
                self.obs = obs
                self.fluxm = flux
                self.insert = insert

class RWRate(RWReweight):
    subcommand = 'rate'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'reweight.h5'
    default_output_file = 'reweight.h5'
    description = '''\ Blah!'''
    def w_postanalysis_stateprobs():
        pass

    def w_postanalysis_reweight(self):
        self.generate_reweight_data()
        pi = self.progress.indicator

        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

        start_pts = range(start_iter, stop_iter, step_iter)


        indices = np.array(range(start_iter-1, stop_iter-1), dtype=np.uint16)

        submit_kwargs = dict(pi=pi, nstates=self.nstates, start_iter=self.start_iter, stop_iter=self.stop_iter, 
                             step_iter=self.step_iter)

        # Due to the way our estimator is written, we're using the same dataset every time.  We're just returning different values.
        submit_kwargs['dataset'] = {'indices': indices}
        submit_kwargs['estimator_kwargs'] = {}
        submit_kwargs['estimator_kwargs'].update(  dict(rows=self.rows,
                                                        cols=self.cols,
                                                        obs=self.obs,
                                                        flux=self.fluxm,
                                                        insert=self.insert,
                                                        bin_last_state_map=np.tile(np.arange(self.nstates, dtype=np.int), self.nbins), 
                                                        bin_state_map=np.repeat(self.state_map[:-1], self.nstates),
                                                        nfbins=self.nfbins,
                                                        state_labels=self.state_labels,
                                                        return_flux=False,
                                                        return_states=False,
                                                        return_color=False,
                                                        state_map=self.state_map,
                                                        nbins=self.nbins))


        # Calculate averages for the simulation, then report, if necessary.
        with pi:

            # The dataset options are what we pass on to the estimator...
            submit_kwargs['dataset'].update(dict(pre_calculated=self.rates))
            avg_rates = self.run_calculation(eval_block=_2D_eval_block, name='Rate Evolution', dim=2, do_averages=True, **submit_kwargs)
            self.output_file.replace_dataset('avg_rates', data=avg_rates[1])

            submit_kwargs['dataset'].update(dict(pre_calculated=self.flux))
            submit_kwargs['estimator_kwargs']['return_flux'] = True
            avg_conditional_fluxes = self.run_calculation(eval_block=_2D_eval_block, name='Flux Evolution', dim=2, do_averages=True, **submit_kwargs)
            self.output_file.replace_dataset('avg_conditional_fluxes', data=avg_conditional_fluxes[1])

        # Now, print them!
        pi.clear()

        # We've returned an average, but it still exists in a timeslice format.  So we need to return the 'last' value.
        self.print_averages(avg_conditional_fluxes[1], '\nfluxes from state to state:', dim=2)
        self.print_averages(avg_rates[1], '\nrates from state to state:', dim=2)

        # Do a bootstrap evolution.
        pi.clear()
        with pi:

            submit_kwargs['dataset'].update(dict(pre_calculated=self.rates))
            submit_kwargs['estimator_kwargs']['return_flux'] = False
            rate_evol = self.run_calculation(eval_block=_2D_eval_block, name='Rate Evolution', dim=2, **submit_kwargs)
            rate_evol['expected'] *= (self.npts - 1)
            self.output_file.replace_dataset('rate_evolution', data=rate_evol, shuffle=True, compression=9)
            pi.clear()

            submit_kwargs['dataset'].update(dict(pre_calculated=self.flux))
            submit_kwargs['estimator_kwargs']['return_flux'] = True
            flux_evol = self.run_calculation(eval_block=_2D_eval_block, name='Conditional Flux Evolution', dim=2, **submit_kwargs)
            flux_evol['expected'] *= (self.npts - 1)
            self.output_file.replace_dataset('conditional_flux_evolution', data=rate_evol, shuffle=True, compression=9)
            pi.clear()

    def go(self):
        self.w_postanalysis_reweight()


class RWStateProbs(RWReweight):
    subcommand = 'stateprobs'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'reweight.h5'
    description = '''\ Blah!'''
    def w_postanalysis_stateprobs(self):
        self.generate_reweight_data()
        pi = self.progress.indicator

        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

        start_pts = range(start_iter, stop_iter, step_iter)


        indices = np.array(range(start_iter-1, stop_iter-1), dtype=np.uint16)

        submit_kwargs = dict(pi=pi, nstates=self.nstates, start_iter=self.start_iter, stop_iter=self.stop_iter, 
                             step_iter=self.step_iter)

        # Due to the way our estimator is written, we're using the same dataset every time.  We're just returning different values.
        submit_kwargs['dataset'] = {'indices': indices}
        submit_kwargs['estimator_kwargs'] = {}
        submit_kwargs['estimator_kwargs'].update(  dict(rows=self.rows,
                                                        cols=self.cols,
                                                        obs=self.obs,
                                                        flux=self.fluxm,
                                                        insert=self.insert,
                                                        bin_last_state_map=np.tile(np.arange(self.nstates, dtype=np.int), self.nbins), 
                                                        bin_state_map=np.repeat(self.state_map[:-1], self.nstates),
                                                        nfbins=self.nfbins,
                                                        state_labels=self.state_labels,
                                                        return_flux=False,
                                                        return_states=False,
                                                        return_color=True,
                                                        state_map=self.state_map,
                                                        nbins=self.nbins))


        # Calculate averages for the simulation, then report, if necessary.
        with pi:

            # The dataset options are what we pass on to the estimator...
            submit_kwargs['dataset'].update(dict(pre_calculated=self.color_prob))
            avg_color_probs = self.run_calculation(eval_block=_1D_eval_block, name='Color (Ensemble) Probability Evolution', dim=1, do_averages=True, **submit_kwargs)
            self.output_file.replace_dataset('avg_color_probs', data=avg_color_probs[1])

            submit_kwargs['dataset'].update(dict(pre_calculated=self.state_prob))
            submit_kwargs['estimator_kwargs']['return_color'] = False
            submit_kwargs['estimator_kwargs']['return_states'] = True
            avg_state_probs = self.run_calculation(eval_block=_1D_eval_block, name='State Probability Evolution', dim=1, do_averages=True, **submit_kwargs)
            self.output_file.replace_dataset('avg_state_probs', data=avg_state_probs[1])

        # Now, print them!
        pi.clear()

        # We've returned an average, but it still exists in a timeslice format.  So we need to return the 'last' value.
        self.print_averages(avg_color_probs[1], '\naverage color probabilities:', dim=1)
        self.print_averages(avg_state_probs[1], '\naverage state probabilities:', dim=1)

        # Do a bootstrap evolution.
        pi.clear()
        with pi:

            submit_kwargs['dataset'].update(dict(pre_calculated=self.color_prob))
            submit_kwargs['estimator_kwargs']['return_color'] = True
            submit_kwargs['estimator_kwargs']['return_states'] = False
            color_evol = self.run_calculation(eval_block=_1D_eval_block, name='Color (Ensemble) Probability Evolution', dim=1, **submit_kwargs)
            #color_evol['expected'] *= (self.npts - 1)
            self.output_file.replace_dataset('color_prob_evolution', data=color_evol, shuffle=True, compression=9)
            pi.clear()

            submit_kwargs['dataset'].update(dict(pre_calculated=self.state_prob))
            submit_kwargs['estimator_kwargs']['return_color'] = False
            submit_kwargs['estimator_kwargs']['return_states'] = True
            state_evol = self.run_calculation(eval_block=_1D_eval_block, name='State Probability Evolution', dim=1, **submit_kwargs)
            #state_evol['expected'] *= (self.npts - 1)
            self.output_file.replace_dataset('state_pop_evolution', data=state_evol, shuffle=True, compression=9)
            pi.clear()

    def go(self):
        self.w_postanalysis_stateprobs()

class RWAll(RWMatrix, RWStateProbs, RWRate):
    subcommand = 'all'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'reweight.h5'

    def go(self):
        # One minor issue; as this stands now, since it's inheriting from all the other classes, it needs
        # a kinetics file to instantiate the other attributes.  We'll need to modify how the loading works, there.
        self.w_postanalysis_matrix()
        self.w_postanalysis_reweight()
        self.w_postanalysis_stateprobs()

# Just a convenience class to average the observables.
class RWAverage(RWStateProbs, RWRate):
    subcommand = 'average'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'reweight.h5'
    default_output_file = 'reweight.h5'

    def go(self):
        self.w_postanalysis_reweight()
        self.w_postanalysis_stateprobs()

class WReweight(WESTMasterCommand, WESTParallelTool):
    prog='w_reweight'
    #subcommands = [AvgTraceSubcommand,AvgMatrixSubcommand]
    subcommands = [RWRate, RWStateProbs, RWAll, RWAverage, RWMatrix]
    subparsers_title = 'direct kinetics analysis schemes'

if __name__ == '__main__':
    WReweight().main()
