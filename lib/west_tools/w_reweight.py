# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
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

import logging

# Let's suppress those numpy warnings.
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

import numpy as np
import scipy.sparse as sp
import h5py

import westpa
from west.data_manager import weight_dtype, n_iter_dtype
#from westtools import (WESTTool, WESTParallelTool, WESTDataReader, IterRangeSelection,
from westtools import (WESTMasterCommand, WESTParallelTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)

from westtools.kinetics_tool import WESTKineticsBase, AverageCommands, generate_future
from westpa import h5io
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

log = logging.getLogger('westtools.w_reweight')

import mclib

from mclib import mcbs_correltime, mcbs_ci_correl

# From postanalysis matrix
from westpa.binning import index_dtype
from westpa.reweight import stats_process, reweight_for_c, FluxMatrix

def _2D_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, mcbs_enable, estimator_kwargs):
    # As our reweighting estimator is a weird function, we can't use the general mclib block.
    results = []
    for istate in range(nstates):
        for jstate in range(nstates):
            if istate == jstate: continue
            estimator_kwargs.update(dict(istate=istate, jstate=jstate, nstates=nstates))

            dataset = { 'indices' : np.array(list(range(start-1, stop-1)), dtype=np.uint16) }
            
            ci_res = mcbs_ci_correl(dataset,estimator=reweight_for_c,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=(lambda x: x[0]), do_correl=do_correl, mcbs_enable=mcbs_enable, estimator_kwargs=estimator_kwargs)
            results.append((name, iblock, istate, jstate, (start,stop) + ci_res))

    return results

def _1D_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, mcbs_enable, estimator_kwargs):
    # As our reweighting estimator is a weird function, we can't use the general mclib block.
    results = []
    for istate in range(nstates):
        # A little hack to make our estimator play nice, as jstate must be there.
        # For 1D datasets (state probabilities, etc), the argument isn't used in our estimator,
        # and so any variable which has the proper type is fine.
        estimator_kwargs.update(dict(istate=istate, jstate=istate, nstates=nstates))

        dataset = { 'indices' : np.array(list(range(start-1, stop-1)), dtype=np.uint16) }
        
        ci_res = mcbs_ci_correl(dataset,estimator=reweight_for_c,
                                alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                subsample=(lambda x: x[0]), do_correl=do_correl, mcbs_enable=mcbs_enable, estimator_kwargs=estimator_kwargs)
        results.append((name, iblock, istate, (start,stop) + ci_res))

    return results

def _pop_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, do_correl, mcbs_enable, estimator_kwargs, **kwargs):
    # As our reweighting estimator is a weird function, we can't use the general mclib block.
    results = []
    # A little hack to make our estimator play nice, as jstate must be there.
    # For 1D datasets (state probabilities, etc), the argument isn't used in our estimator,
    # and so any variable which has the proper type is fine.
    estimator_kwargs.update(dict(istate=0, jstate=0, nstates=nstates))

    estimator_kwargs.update({ 'indices' : np.array(list(range(start-1, stop-1)), dtype=np.uint16), 'stride': 1  })
    #cpdef reweight_for_c(rows, cols, obs, flux, insert, indices, nstates, nbins, state_labels, state_map, nfbins, istate, jstate, stride, bin_last_state_map, bin_state_map, return_obs, obs_threshold=1):
    #['nbins', 'rows', 'state_map', 'nfbins', 'bin_state_map', 'cols', 'jstate', 'flux', 'istate', 'return_obs', 'bin_last_state_map', 'insert', 'nstates', 'indices', 'state_labels', 'obs']
    #cpdef reweight_for_c(stride, obs_threshold=1):
    
    ci_res = reweight_for_c(**estimator_kwargs)[...].reshape(-1,nstates).sum(axis=1)

    return ci_res

class RWMatrix(WESTKineticsBase, FluxMatrix):
    subcommand = 'init'
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

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "reweight.h5") contains the
following datasets:

  ``/bin_populations`` [window, bin]
     The reweighted populations of each bin based on windows. Bins contain
     one color each, so to recover the original un-colored spatial bins,
     one must sum over all states.

  ``/iterations`` [iteration]
    *(Structured -- see below)*  Sparse matrix data from each
    iteration.  They are reconstructed and averaged within the
    w_reweight {kinetics/probs} routines so that observables may
    be calculated.  Each group contains 4 vectors of data:
      
      flux
        *(Floating-point)* The weight of a series of flux events
      cols
        *(Integer)* The bin from which a flux event began.
      cols
        *(Integer)* The bin into which the walker fluxed.
      obs
        *(Integer)* How many flux events were observed during this
        iteration.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''
    
    def __init__(self, parent):
        super(RWMatrix, self).__init__(parent)
        self.parent = parent
        self.assignments_file = None
        

    def add_args(self, parser):
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

                                       
    def process_args(self, args):
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        # Force a build of the transition matrix at the iteration level.
        self.sampling_frequency = 'iteration' if self.assignments_file.attrs['subsampled'] == True else args.sampling_frequency

    def go(self):
        # This function exists in FluxMatrix, and is not currently portable.
        # That is, it assumes all properties are properly initialized (self.assignments_file, etc)
        # TO DO : make it portable.
        pi = self.progress.indicator
        with pi:
            self.w_postanalysis_matrix()

class RWReweight(AverageCommands):
    help_text = 'Parent class for all reweighting routines, as they all use the same estimator code.'

    def __init__(self, parent):
        super(RWReweight, self).__init__(parent)
        self.parent = parent
        self.assignments_file = None

    def add_args(self, parser):
        cogroup = parser.add_argument_group('calculation options')
        cogroup.add_argument('--obs-threshold', type=int, default=1,
                             help='''The minimum number of observed transitions between two states i and j necessary to include
                             fluxes in the reweighting estimate''')

    def process_args(self, args):
        self.obs_threshold = args.obs_threshold

    def accumulate_statistics(self,start_iter,stop_iter):
        ''' 
        This function pulls previously generated flux matrix data into memory.
        The data is assumed to exist within an HDF5 file that is available
        as a property.
        The data is kept as a single dimensional numpy array to use with the
        cython estimator.
        '''
        # This is designed to pull in the flux data...
        rows = []
        cols = []
        obs = []
        flux = []
        insert = [0]*(start_iter)

        # Actually, I'm not sure we need to start this at start_iter...
        # ... as it's keyed to the iteration, we need to make sure that the index
        # matches with the iteration (particularly for 'insert').
        # It's just easier to load all the data, although we could just start insert as a list of length
        # start_iter.
        for iiter in range(start_iter, stop_iter):
            iter_grp = self.kinetics_file['iterations']['iter_{:08d}'.format(iiter)]

            rows.append(iter_grp['rows'][...])
            cols.append(iter_grp['cols'][...])
            obs.append(iter_grp['obs'][...])
            flux.append(iter_grp['flux'][...])
            # 'insert' is the insertion point for each iteration; that is,
            # at what point do we look into the list for iteration X?
            insert.append(iter_grp['rows'][...].shape[0] + insert[-1])
        self.rows = np.concatenate(rows)
        self.cols = np.concatenate(cols)
        self.obs = np.concatenate(obs)
        self.flux = np.concatenate(flux)
        assert insert[-1] == len(self.rows)
        self.insert = np.array(insert, dtype=np.intc)

    def generate_reweight_data(self):
        ''' 
        This function ensures all the appropriate files are loaded, sets
        appropriate attributes necessary for all calling functions/children,
        and then calls the function to load in the flux matrix data.
        '''

        self.open_files()
        self.open_assignments()

        self.nfbins = self.kinetics_file.attrs['nrows']
        self.npts = self.kinetics_file.attrs['npts']
        self.state_map = self.assignments_file['state_map'][...]
        self.state_labels = self.assignments_file['state_labels'][...]

        # Copying this over to make it more convenient...
        self.output_file.replace_dataset('state_labels', data=self.assignments_file['state_labels'][...])
        # ... as we'll need it for ploterr.

        assert self.nstates == len(self.state_labels)
        assert self.nfbins == self.nbins * self.nstates

        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

        # This function call loads up all the flux matrix data and stores it in memory.
        # We just give our cython estimator all the information it needs beforehand, as
        # the nature of the bootstrap means we can't know what iterations we'll actually need before
        # the routine is called.

        self.accumulate_statistics(start_iter,stop_iter)

class RWRate(RWReweight):
    subcommand = 'kinetics'
    help_text = 'Generates rate and flux values from a WESTPA simulation via reweighting.'
    default_kinetics_file = 'reweight.h5'
    default_output_file = 'reweight.h5'
    description = '''\
Calculate average rates from weighted ensemble data using the postanalysis
reweighting scheme. Bin assignments (usually "assign.h5") and pre-calculated 
iteration flux matrices (usually "reweight.h5") data files must have been 
previously generated using w_reweight matrix (see "w_assign --help" and 
"w_reweight init --help" for information on generating these files).

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------
The output file (-o/--output, usually "kinrw.h5") contains the following
dataset:

  /avg_rates [state,state]
    (Structured -- see below) State-to-state rates based on entire window of
    iterations selected.

  /avg_total_fluxes [state]
    (Structured -- see below) Total fluxes into each state based on entire
    window of iterations selected.
    
  /avg_conditional_fluxes [state,state]
    (Structured -- see below) State-to-state fluxes based on entire window of
    iterations selected.

If --evolution-mode is specified, then the following additional datasets are
available:

  /rate_evolution [window][state][state]
    (Structured -- see below). State-to-state rates based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive 
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.
    
  /target_flux_evolution [window,state]
    (Structured -- see below). Total flux into a given macro state based on
    windows of iterations of varying width, as in /rate_evolution.
    
  /conditional_flux_evolution [window,state,state]
    (Structured -- see below). State-to-state fluxes based on windows of
    varying width, as in /rate_evolution.
    
The structure of these datasets is as follows:

  iter_start
    (Integer) Iteration at which the averaging window begins (inclusive).
    
  iter_stop
    (Integer) Iteration at which the averaging window ends (exclusive).
    
  expected
    (Floating-point) Expected (mean) value of the observable as evaluated within
    this window, in units of inverse tau.
    
  ci_lbound
    (Floating-point) Lower bound of the confidence interval of the observable
    within this window, in units of inverse tau.
    
  ci_ubound
    (Floating-point) Upper bound of the confidence interval of the observable
    within this window, in units of inverse tau.

  stderr
    (Floating-point) The standard error of the mean of the observable
    within this window, in units of inverse tau.
    
  corr_len
    (Integer) Correlation length of the observable within this window, in units
    of tau.

Each of these datasets is also stamped with a number of attributes:

  mcbs_alpha
    (Floating-point) Alpha value of confidence intervals. (For example, 
    *alpha=0.05* corresponds to a 95% confidence interval.)

  mcbs_nsets
    (Integer) Number of bootstrap data sets used in generating confidence
    intervals.
    
  mcbs_acalpha
    (Floating-point) Alpha value for determining correlation lengths.
   

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
    '''

    def __init__(self, parent):
        super(RWRate, self).__init__(parent)
        self.parent = parent
        self.assignments_file = None

    def w_postanalysis_reweight(self):
        ''' 
        This function ensures the data is ready to send in to the estimator and
        the bootstrapping routine, then does so. Much of this is simply setting 
        up appropriate args and kwargs, then passing them in to the 
        'run_calculation' function, which sets up future objects to send to the
        work manager. The results are returned, and then written to the 
        appropriate HDF5 dataset.
        This function is specific for the rates and fluxes from the reweighting 
        method.
        '''
        self.generate_reweight_data()
        pi = self.progress.indicator

        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step
        start_pts = list(range(start_iter, stop_iter, step_iter))
        indices = np.array(list(range(start_iter-1, stop_iter-1)), dtype=np.uint16)
        submit_kwargs = dict(pi=pi, nstates=self.nstates, start_iter=self.start_iter, stop_iter=self.stop_iter, 
                             step_iter=self.step_iter)

        # Due to the way our estimator is written, we're using the same dataset every time.  We're just returning different values.
        submit_kwargs['dataset'] = {'indices': indices}
        submit_kwargs['estimator_kwargs'] = {}
        submit_kwargs['estimator_kwargs'].update(  dict(rows=self.rows,
                                                        cols=self.cols,
                                                        obs=self.obs,
                                                        flux=self.flux,
                                                        insert=self.insert,
                                                        bin_last_state_map=np.tile(np.arange(self.nstates, dtype=np.int), self.nbins), 
                                                        bin_state_map=np.repeat(self.state_map[:-1], self.nstates),
                                                        nfbins=self.nfbins,
                                                        state_labels=self.state_labels,
                                                        return_obs='R', # Set to a default, here, but we explicitly set it later.
                                                        state_map=self.state_map,
                                                        nbins=self.nbins))


        # Calculate averages for the simulation, then report, if necessary.

        # The dataset options are what we pass on to the estimator...
        submit_kwargs['estimator_kwargs']['return_obs'] = 'R'
        avg_rates = self.run_calculation(eval_block=_2D_eval_block, name='Average Rates', dim=2, do_averages=True, **submit_kwargs)
        avg_rates['expected'] *= (self.npts - 1)
        avg_rates['ci_ubound'] *= (self.npts - 1)
        avg_rates['ci_lbound'] *= (self.npts - 1)
        self.output_file.replace_dataset('avg_rates', data=avg_rates[1])

        submit_kwargs['estimator_kwargs']['return_obs'] = 'F'
        avg_conditional_fluxes = self.run_calculation(eval_block=_2D_eval_block, name='Average Flux', dim=2, do_averages=True, **submit_kwargs)
        avg_conditional_fluxes['expected'] *= (self.npts - 1)
        avg_conditional_fluxes['ci_ubound'] *= (self.npts - 1)
        avg_conditional_fluxes['ci_lbound'] *= (self.npts - 1)
        self.output_file.replace_dataset('avg_conditional_fluxes', data=avg_conditional_fluxes[1])

        # Now, print them!

        # We've returned an average, but it still exists in a timeslice format.  So we need to return the 'last' value.
        if self.display_averages:
            self.print_averages(avg_conditional_fluxes[1], '\nfluxes from state to state:', dim=2)
            self.print_averages(avg_rates[1], '\nrates from state to state:', dim=2)

        # Do a bootstrap evolution.

        submit_kwargs['estimator_kwargs']['return_obs'] = 'R'
        rate_evol = self.run_calculation(eval_block=_2D_eval_block, name='Rate Evolution', dim=2, **submit_kwargs)
        rate_evol['expected'] *= (self.npts - 1)
        rate_evol['ci_ubound'] *= (self.npts - 1)
        rate_evol['ci_lbound'] *= (self.npts - 1)
        self.output_file.replace_dataset('rate_evolution', data=rate_evol, shuffle=True, compression=9)

        submit_kwargs['estimator_kwargs']['return_obs'] = 'F'
        flux_evol = self.run_calculation(eval_block=_2D_eval_block, name='Conditional Flux Evolution', dim=2, **submit_kwargs)
        flux_evol['expected'] *= (self.npts - 1)
        flux_evol['ci_ubound'] *= (self.npts - 1)
        flux_evol['ci_lbound'] *= (self.npts - 1)
        self.output_file.replace_dataset('conditional_flux_evolution', data=rate_evol, shuffle=True, compression=9)


    def go(self):
        pi = self.progress.indicator
        with pi:
            self.w_postanalysis_reweight()


class RWStateProbs(RWReweight):
    subcommand = 'probs'
    help_text = 'Calculates color and state probabilities via reweighting.'
    default_kinetics_file = 'reweight.h5'
    description = '''\
Calculate average populations from weighted ensemble data using the postanalysis
reweighting scheme. Bin assignments (usually "assign.h5") and pre-calculated 
iteration flux matrices (usually "reweight.h5") data files must have been 
previously generated using w_reweight matrix (see "w_assign --help" and 
"w_reweight init --help" for information on generating these files).

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, usually "direct.h5") contains the following
dataset:

  /avg_state_probs [state]
    (Structured -- see below) Population of each state across entire
    range specified.

  /avg_color_probs [state]
    (Structured -- see below) Population of each ensemble across entire
    range specified.

If --evolution-mode is specified, then the following additional datasets are
available:

  /state_pop_evolution [window][state]
    (Structured -- see below). State populations based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive 
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.

  /color_prob_evolution [window][state]
    (Structured -- see below). Ensemble populations based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive 
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.
    
The structure of these datasets is as follows:

  iter_start
    (Integer) Iteration at which the averaging window begins (inclusive).
    
  iter_stop
    (Integer) Iteration at which the averaging window ends (exclusive).
    
  expected
    (Floating-point) Expected (mean) value of the observable as evaluated within
    this window, in units of inverse tau.
    
  ci_lbound
    (Floating-point) Lower bound of the confidence interval of the observable
    within this window, in units of inverse tau.
    
  ci_ubound
    (Floating-point) Upper bound of the confidence interval of the observable
    within this window, in units of inverse tau.

  stderr
    (Floating-point) The standard error of the mean of the observable
    within this window, in units of inverse tau.
    
  corr_len
    (Integer) Correlation length of the observable within this window, in units
    of tau.


Each of these datasets is also stamped with a number of attributes:

  mcbs_alpha
    (Floating-point) Alpha value of confidence intervals. (For example, 
    *alpha=0.05* corresponds to a 95% confidence interval.)

  mcbs_nsets
    (Integer) Number of bootstrap data sets used in generating confidence
    intervals.
    
  mcbs_acalpha
    (Floating-point) Alpha value for determining correlation lengths.
   

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''    
    def w_postanalysis_stateprobs(self):
        ''' 
        This function ensures the data is ready to send in to the estimator and
        the bootstrapping routine, then does so. Much of this is simply setting 
        up appropriate args and kwargs, then passing them in to the 
        'run_calculation' function, which sets up future objects to send to the
        work manager. The results are returned, and then written to the 
        appropriate HDF5 dataset.
        This function is specific for the color (steady-state) and macrostate 
        probabilities from the reweighting method.
        '''
        self.generate_reweight_data()
        pi = self.progress.indicator

        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step
        start_pts = list(range(start_iter, stop_iter, step_iter))
        indices = np.array(list(range(start_iter-1, stop_iter-1, step_iter)), dtype=np.uint16)
        submit_kwargs = dict(pi=pi, nstates=self.nstates, start_iter=self.start_iter, stop_iter=self.stop_iter, 
                             step_iter=self.step_iter)

        # Due to the way our estimator is written, we're using the same dataset every time.  We're just returning different values.
        submit_kwargs['dataset'] = {'indices': indices}
        submit_kwargs['estimator_kwargs'] = {}
        submit_kwargs['estimator_kwargs'].update(  dict(rows=self.rows,
                                                        cols=self.cols,
                                                        obs=self.obs,
                                                        flux=self.flux,
                                                        insert=self.insert,
                                                        bin_last_state_map=np.tile(np.arange(self.nstates, dtype=np.int), self.nbins), 
                                                        bin_state_map=np.repeat(self.state_map[:-1], self.nstates),
                                                        nfbins=self.nfbins,
                                                        state_labels=self.state_labels,
                                                        return_obs='C', # Set to a default, but we explicitly set it later.
                                                        state_map=self.state_map,
                                                        nbins=self.nbins))


        # Calculate averages for the simulation, then report, if necessary.

        # The dataset options are what we pass on to the estimator...
        submit_kwargs['estimator_kwargs']['return_obs'] = 'C'
        avg_color_probs = self.run_calculation(eval_block=_1D_eval_block, name='Average Color (Ensemble) Probability', dim=1, do_averages=True, **submit_kwargs)
        self.output_file.replace_dataset('avg_color_probs', data=avg_color_probs[1])

        submit_kwargs['estimator_kwargs']['return_obs'] = 'S'
        avg_state_probs = self.run_calculation(eval_block=_1D_eval_block, name='Average State Probability', dim=1, do_averages=True, **submit_kwargs)
        self.output_file.replace_dataset('avg_state_probs', data=avg_state_probs[1])

        # Now, print them!

        # We've returned an average, but it still exists in a timeslice format.  So we need to return the 'last' value.
        if self.display_averages:
            self.print_averages(avg_color_probs[1], '\naverage color probabilities:', dim=1)
            self.print_averages(avg_state_probs[1], '\naverage state probabilities:', dim=1)

        # Do a bootstrap evolution.

        submit_kwargs['estimator_kwargs']['return_obs'] = 'C'
        color_evol = self.run_calculation(eval_block=_1D_eval_block, name='Color (Ensemble) Probability Evolution', dim=1, **submit_kwargs)
        self.output_file.replace_dataset('color_prob_evolution', data=color_evol, shuffle=True, compression=9)

        submit_kwargs['estimator_kwargs']['return_obs'] = 'S'
        state_evol = self.run_calculation(eval_block=_1D_eval_block, name='State Probability Evolution', dim=1, **submit_kwargs)
        self.output_file.replace_dataset('state_pop_evolution', data=state_evol, shuffle=True, compression=9)

        # Little different, here.  Do a non-bootstrap evolution.
        # Just use the work manager and the estimator without the bootstrap code, which is problematic for a dataset like this.
        # We'll be adding this in later.
        if False:
            submit_kwargs['estimator_kwargs']['return_obs'] = 'P'
            submit_kwargs['pi'] = None
            futures = []
            bin_pops = []
            for iblock, start in enumerate(start_pts):
                stop = min(start+step_iter, stop_iter)
                if self.evolution_mode == 'cumulative' or do_averages == True:
                    windowsize = int(self.evol_window_frac * (stop - start_iter))
                    block_start = max(start_iter, stop - windowsize)
                else: # self.evolution_mode == 'blocked'
                    block_start = start
                future_kwargs = dict(iblock=iblock, start=block_start, stop=stop,
                                     #nstates=self.nstates,
                                     mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                     mcbs_acalpha=self.mcbs_acalpha,
                                     do_correl=self.do_correl,name='Bin Population Evolution',
                                     mcbs_enable=self.mcbs_enable,
                                     data_input={},
                                     **submit_kwargs)
                #print(future_kwargs)
                futures.append(generate_future(self.work_manager, 'Bin Pop Evolution', _pop_eval_block, future_kwargs))
            #bin_evol = self.run_calculation(eval_block=_pop_eval_block, name='Bin Population Evolution', dim=1, **submit_kwargs)
            #bin_evol = bin_evol.reshape(-1, 2)
            for future in self.work_manager.as_completed(futures):
                #pi.progress += iblock / step_iter
                bin_pops.append(future.get_result(discard=True))
            hist = self.output_file.replace_dataset('histograms', data=bin_pops, shuffle=True, compression=9)
            hist.attrs['iter_start'] = start_iter
            hist.attrs['iter_stop'] = stop_iter
            import ast
            binbounds = []
            midpoints = []
            for i in self.assignments_file['bin_labels']:
                binbounds.append(ast.literal_eval(i)[0][0])
            # This is probably crap, so...
            binbounds.append(ast.literal_eval(self.assignments_file['bin_labels'][-1])[0][1])
            #binbounds.append(binbounds[-1]+1)
            for i in range(1, len(binbounds)):
                midpoints.append(((binbounds[i] - binbounds[i-1])/2)+binbounds[i-1])
            self.output_file.replace_dataset('binbounds_0', data=binbounds, shuffle=True, compression=9)
            self.output_file.replace_dataset('midpoints_0', data=midpoints, shuffle=True, compression=9)
            self.output_file.replace_dataset('n_iter', data=list(range(start_iter,stop_iter)), shuffle=True, compression=9)

    def go(self):
        pi = self.progress.indicator
        with pi:
            self.w_postanalysis_stateprobs()

class RWAll(RWMatrix, RWStateProbs, RWRate):
    subcommand = 'all'
    help_text = 'Runs the full suite, including the generation of the flux matrices.'
    default_kinetics_file = 'reweight.h5'
    default_output_file = 'reweight.h5'
    description = '''\
A convenience function to run init/kinetics/probs. Bin assignments, 
including macrostate definitions, are required. (See 
"w_assign --help" for more information).

For more information on the individual subcommands this subs in for, run
w_reweight {init/kinetics/probs} --help.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''    

    def go(self):
        pi = self.progress.indicator
        with pi:
            self.w_postanalysis_matrix()
            self.w_postanalysis_reweight()
            self.w_postanalysis_stateprobs()

# Just a convenience class to average the observables.
class RWAverage(RWStateProbs, RWRate):
    subcommand = 'average'
    help_text = 'Averages and returns fluxes, rates, and color/state populations.'
    default_kinetics_file = 'reweight.h5'
    default_output_file = 'reweight.h5'
    description = '''\
A convenience function to run kinetics/probs. Bin assignments, 
including macrostate definitions, are required. (See 
"w_assign --help" for more information).

For more information on the individual subcommands this subs in for, run
w_reweight {kinetics/probs} --help.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''    

    def go(self):
        pi = self.progress.indicator
        with pi:
            self.w_postanalysis_reweight()
            self.w_postanalysis_stateprobs()

class WReweight(WESTMasterCommand, WESTParallelTool):
    prog='w_reweight'
    subcommands = [RWMatrix, RWAverage, RWRate, RWStateProbs, RWAll]
    subparsers_title = 'reweighting kinetics analysis scheme'

if __name__ == '__main__':
    WReweight().main()
