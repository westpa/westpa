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

import sys, random, math
import numpy, h5py
from h5py import h5s

import westpa
from west.data_manager import weight_dtype, n_iter_dtype
from westtools import (WESTMasterCommand, WESTParallelTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)
from westpa import h5io
from westpa.kinetics import labeled_flux_to_rate, sequence_macro_flux_to_rate, sequence_macro_flux_to_rate_bs
from westpa.kinetics.matrates import get_macrostate_rates

import mclib
from mclib import mcbs_correltime, mcbs_ci_correl, mcbs_ci_correl_rw


log = logging.getLogger('westtools.w_kinavg')

from westtools.dtypes import iter_block_ci_dtype as ci_dtype

class KinAvgSubcommands(WESTSubcommand):
    '''Common argument processing for w_kinavg subcommands'''
    
    def __init__(self, parent):
        super(KinAvgSubcommands,self).__init__(parent)
        
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
        
        # self.default_kinetics_file will be picked up as a class attribute from the appropriate subclass        
        iogroup.add_argument('-k', '--kinetics', default=self.default_kinetics_file,
                            help='''Populations and transition rates are stored in KINETICS
                            (default: %(default)s).''')
        iogroup.add_argument('-o', '--output', dest='output', default='kinavg.h5',
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
        cgroup.add_argument('--nsets', type=int,
                             help='''Use NSETS samples for bootstrapping (default: chosen based on ALPHA)''')
        
        cogroup = parser.add_argument_group('calculation options')
        cogroup.add_argument('-e', '--evolution-mode', choices=['cumulative', 'blocked', 'none'], default='none',
                             help='''How to calculate time evolution of rate estimates.
                             ``cumulative`` evaluates rates over windows starting with --start-iter and getting progressively
                             wider to --stop-iter by steps of --step-iter.
                             ``blocked`` evaluates rates over windows of width --step-iter, the first of which begins at
                             --start-iter.
                             ``none`` (the default) disables calculation of the time evolution of rate estimates.''')
        cogroup.add_argument('--window-frac', type=float, default=1.0,
                             help='''Fraction of iterations to use in each window when running in ``cumulative`` mode.
                             The (1 - frac) fraction of iterations will be discarded from the start of each window.''')

        mgroup = parser.add_argument_group('misc options')
        mgroup.add_argument('--disable-averages', '-da', dest='display_averages', action='store_false',
                             help='''Whether or not the averages should be printed to the console (set to FALSE if flag is used).''')
        
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
        self.correl = args.correl
        self.mcbs_alpha = args.alpha
        self.mcbs_acalpha = args.acalpha if args.acalpha else self.mcbs_alpha
        self.mcbs_nsets = args.nsets if args.nsets else mclib.get_bssize(self.mcbs_alpha)

        self.display_averages = args.display_averages
        
        self.evolution_mode = args.evolution_mode
        self.evol_window_frac = args.window_frac
        if self.evol_window_frac <= 0 or self.evol_window_frac > 1:
            raise ValueError('Parameter error -- fractional window defined by --window-frac must be in (0,1]')

def _eval_block(iblock, start, stop, nstates, total_fluxes, cond_fluxes, pops, rates, mcbs_alpha, mcbs_nsets, mcbs_acalpha, correl):
    results = [[],[],[]]
    # results are target fluxes, conditional fluxes, rates
    for istate in xrange(nstates):
        kwargs = { 'istate' : istate }
        dataset = {'dataset': total_fluxes[:, istate]}
        ci_res = mcbs_ci_correl_rw(dataset,estimator=(lambda stride, dataset: numpy.mean(dataset)),
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=numpy.mean, pre_calculated=dataset['dataset'], correl=correl)
        results[0].append((iblock,istate,(start,stop)+ci_res))
        
        for jstate in xrange(nstates):
            if istate == jstate: continue
            kwargs = { 'istate' : istate, 'jstate': jstate }
            dataset = {'dataset': cond_fluxes[:, istate, jstate]}
            ci_res = mcbs_ci_correl_rw(dataset,estimator=(lambda stride, dataset: numpy.mean(dataset)),
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=numpy.mean, pre_calculated=dataset['dataset'], correl=correl)
            results[1].append((iblock, istate, jstate, (start,stop) + ci_res))
            
            # macro_flux_to_rate_bs needs the following:
            # dataset, pops, istate, jstate
            kwargs = { 'istate' : istate, 'jstate': jstate }
            dataset = {'dataset': cond_fluxes[:, istate, jstate], 'pops': pops}
            ci_res = mcbs_ci_correl_rw(dataset,estimator=sequence_macro_flux_to_rate_bs,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=numpy.mean, pre_calculated=rates[:,istate,jstate], correl=correl, **kwargs)
            results[2].append((iblock, istate, jstate, (start,stop) + ci_res))

    return results

class AvgTraceSubcommand(KinAvgSubcommands):
    subcommand = 'trace'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'kintrace.h5'
    
    def __init__(self, parent):
        super(AvgTraceSubcommand,self).__init__(parent)
                        
    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            nstates = self.assignments_file.attrs['nstates']
            nbins = self.assignments_file.attrs['nbins']
            state_labels = self.assignments_file['state_labels'][...]
            assert nstates == len(state_labels)
            start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step
            
            pi.new_operation('Reading data')
            cond_fluxes = h5io.IterBlockedDataset(self.kinetics_file['conditional_fluxes'])
            cond_fluxes.cache_data()
            total_fluxes = h5io.IterBlockedDataset(self.kinetics_file['total_fluxes'])
            pops = h5io.IterBlockedDataset(self.assignments_file['labeled_populations'])
            pops.cache_data()
            pops.data = pops.data.sum(axis=2)
            
            rates = h5io.IterBlockedDataset.empty_like(cond_fluxes)
            rates.data = sequence_macro_flux_to_rate(cond_fluxes.data, pops.data[:nstates,:nbins])
            
            avg_total_fluxes = numpy.zeros((nstates,), dtype=ci_dtype)
            avg_conditional_fluxes = numpy.zeros((nstates,nstates), dtype=ci_dtype)
            avg_rates = numpy.zeros((nstates,nstates), dtype=ci_dtype)
            
            # Calculate overall average rates
            # Replacing the older structure with the futures codepath, because
            # a lot of the way the operations work has changed, and it seems cleaner to use one code path, anyway.
            pi.new_operation('Averaging fluxes/rates', nstates)
            futures = []
            future = self.work_manager.submit(_eval_block, kwargs=dict(iblock=0, start=start_iter, stop=stop_iter,
                                                                       nstates=nstates,
                                                                       total_fluxes=total_fluxes.iter_slice(start_iter,stop_iter),
                                                                       cond_fluxes = cond_fluxes.iter_slice(start_iter,stop_iter),
                                                                       pops=pops.iter_slice(start_iter,stop_iter),
                                                                       rates=rates.iter_slice(start_iter,stop_iter),
                                                                       mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                                                       mcbs_acalpha=self.mcbs_acalpha,
                                                                       correl=self.correl))
            futures.append(future)
            
            for future in self.work_manager.as_completed(futures):
                target_results, condflux_results, rate_results = future.get_result(discard=True)
                for result in target_results:
                    iblock,istate,ci_result = result
                    avg_total_fluxes[istate] = ci_result
                    
                for result in condflux_results:
                    iblock,istate,jstate,ci_result = result
                    avg_conditional_fluxes[istate, jstate] = ci_result
                
                for result in rate_results:
                    iblock, istate, jstate, ci_result = result 
                    avg_rates[istate, jstate] = ci_result
                    
            pi.new_operation('Saving averages')
            self.output_file['avg_rates'] = avg_rates
            self.output_file['avg_conditional_fluxes'] = avg_conditional_fluxes
            self.output_file['avg_total_fluxes'] = avg_total_fluxes
            for ds in ('avg_rates', 'avg_conditional_fluxes', 'avg_total_fluxes'):
                self.stamp_mcbs_info(self.output_file[ds])

            self.output_file['state_labels'] = state_labels
            maxlabellen = max(map(len,state_labels))
            pi.clear()
            
            if self.display_averages:
                print('fluxes into macrostates:')
                for istate in xrange(nstates):
                    print('{:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e}) * tau^-1'
                          .format(state_labels[istate],
                                  avg_total_fluxes['expected'][istate],
                                  avg_total_fluxes['ci_lbound'][istate],
                                  avg_total_fluxes['ci_ubound'][istate],
                                  maxlabellen=maxlabellen))

                print('\nfluxes from state to state:')
                for istate in xrange(nstates):
                    for jstate in xrange(nstates):
                        if istate == jstate: continue
                        print('{:{maxlabellen}s} -> {:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e}) * tau^-1'
                              .format(state_labels[istate], state_labels[jstate],
                                      avg_conditional_fluxes['expected'][istate,jstate],
                                      avg_conditional_fluxes['ci_lbound'][istate,jstate],
                                      avg_conditional_fluxes['ci_ubound'][istate,jstate],
                                      maxlabellen=maxlabellen))
                   
                print('\nrates from state to state:')
                for istate in xrange(nstates):
                    for jstate in xrange(nstates):
                        if istate == jstate: continue
                        print('{:{maxlabellen}s} -> {:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e}) * tau^-1'
                              .format(state_labels[istate], state_labels[jstate],
                                      avg_rates['expected'][istate,jstate],
                                      avg_rates['ci_lbound'][istate,jstate],
                                      avg_rates['ci_ubound'][istate,jstate],
                                      maxlabellen=maxlabellen))
            
            # skip evolution if not requested
            if self.evolution_mode == 'none' or not step_iter: return
            
            start_pts = range(start_iter, stop_iter, step_iter)
            target_evol = numpy.zeros((len(start_pts), nstates), dtype=ci_dtype)
            flux_evol = numpy.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
            rate_evol = numpy.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
            all_items = numpy.arange(1,len(start_pts)+1)
            bootstrap_length = 0.5*(len(start_pts)*(len(start_pts)+1)) - numpy.delete(all_items, numpy.arange(1, len(start_pts)+1, step_iter))
            if self.mcbs_enable == True:
                pi.new_operation('Calculating flux/rate evolution', bootstrap_length[0])
                futures = []
                for iblock, start in enumerate(start_pts):
                    stop = min(start+step_iter, stop_iter)
                    if self.evolution_mode == 'cumulative':
                        windowsize = int(self.evol_window_frac * (stop - start_iter))
                        block_start = max(start_iter, stop - windowsize)
                    else: # self.evolution_mode == 'blocked'
                        block_start = start
                    
                    future = self.work_manager.submit(_eval_block, kwargs=dict(iblock=iblock, start=block_start, stop=stop,
                                                                               nstates=nstates,
                                                                               total_fluxes=total_fluxes.iter_slice(block_start,stop),
                                                                               cond_fluxes = cond_fluxes.iter_slice(block_start,stop),
                                                                               pops=pops.iter_slice(block_start,stop),
                                                                               rates=rates.iter_slice(block_start,stop),
                                                                               mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                                                               mcbs_acalpha=self.mcbs_acalpha,
                                                                               correl=self.correl))
                    futures.append(future)
                
                for future in self.work_manager.as_completed(futures):
                    pi.progress += iblock / step_iter
                    target_results, condflux_results, rate_results = future.get_result(discard=True)
                    for result in target_results:
                        iblock,istate,ci_result = result
                        target_evol[iblock,istate] = ci_result
                        
                    for result in condflux_results:
                        iblock,istate,jstate,ci_result = result
                        flux_evol[iblock,istate, jstate] = ci_result
                    
                    for result in rate_results:
                        iblock, istate, jstate, ci_result = result 
                        rate_evol[iblock, istate, jstate] = ci_result

            else:
                for iblock, start in enumerate(start_pts):
                    stop = min(start+step_iter, stop_iter)
                    if self.evolution_mode == 'cumulative':
                        windowsize = int(self.evol_window_frac * (stop - start_iter))
                        block_start = max(start_iter, stop - windowsize)
                    else: # self.evolution_mode == 'blocked'
                        block_start = start
                    for istate in xrange(nstates):
                        for jstate in xrange(nstates):
                            flux_evol[iblock,istate,jstate]['ci_ubound'] = 0.0
                            flux_evol[iblock,istate,jstate]['ci_lbound'] = 0.0
                            flux_evol[iblock,istate,jstate]['expected'] = numpy.mean(cond_fluxes.iter_slice(block_start, stop)[:, istate, jstate])
                            rate_evol[iblock,istate,jstate]['ci_ubound'] = 0.0
                            rate_evol[iblock,istate,jstate]['ci_lbound'] = 0.0
                            rate_evol[iblock,istate,jstate]['expected'] = rates.iter_slice(block_start, stop)[-1, istate, jstate]

            df_ds = self.output_file.create_dataset('conditional_flux_evolution', data=flux_evol, shuffle=True, compression=9)
            tf_ds = self.output_file.create_dataset('target_flux_evolution', data=target_evol, shuffle=True, compression=9)
            cp_ds = self.output_file.create_dataset('color_prob_evolution', data=pops.data, shuffle=True, compression=9)
            rate_ds = self.output_file.create_dataset('rate_evolution', data=rate_evol, shuffle=True, compression=9)
            
            for ds in (df_ds, tf_ds, rate_ds):
                self.stamp_mcbs_info(ds)

class WKinAvg(WESTMasterCommand, WESTParallelTool):
    prog='w_kinavg'
    #subcommands = [AvgTraceSubcommand,AvgMatrixSubcommand]
    subcommands = [AvgTraceSubcommand]
    subparsers_title = 'kinetics analysis schemes'
    description = '''\
Calculate average rates and associated errors from weighted ensemble data. Bin
assignments (usually "assignments.h5") and kinetics data (usually
"kintrace.h5" or "kinmat.h5") data files must have been previously generated
(see "w_assign --help" and "w_kinetics --help" for information on generating
these files).

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, usually "kinavg.h5") contains the following
dataset:

  /avg_rates [state,state]
    (Structured -- see below) State-to-state rates based on entire window of
    iterations selected.

For trace mode, the following additional datasets are generated:

  /avg_total_fluxes [state]
    (Structured -- see below) Total fluxes into each state based on entire
    window of iterations selected.
    
  /avg_conditional_fluxes [state,state]
    (Structured -- see below) State-to-state fluxes based on entire window of
    iterations selected.

If --evolution-mode is specified, then the following additional dataset is
available:

  /rate_evolution [window][state][state]
    (Structured -- see below). State-to-state rates based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive 
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.
    
If --evolution-mode is specified in trace mode, the following additional
datasets are available:

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
    (Floating-point) Expected (mean) value of the rate as evaluated within
    this window, in units of inverse tau.
    
  ci_lbound
    (Floating-point) Lower bound of the confidence interval on the rate
    within this window, in units of inverse tau.
    
  ci_ubound
    (Floating-point) Upper bound of the confidence interval on the rate 
    within this window, in units of inverse tau.
    
  corr_len
    (Integer) Correlation length of the rate within this window, in units
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

if __name__ == '__main__':
    WKinAvg().main()
