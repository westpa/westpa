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

# Let's suppress those numpy warnings.
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


import sys, random, math
import numpy, h5py
from h5py import h5s

import westpa
from west.data_manager import weight_dtype, n_iter_dtype, seg_id_dtype
from westtools import (WESTMasterCommand, WESTParallelTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)
from westpa import h5io
from westpa.kinetics import labeled_flux_to_rate, sequence_macro_flux_to_rate, sequence_macro_flux_to_rate_bs
from westpa.kinetics.matrates import get_macrostate_rates

import mclib
from mclib import mcbs_correltime, mcbs_ci_correl, mcbs_ci_correl_rw, _1D_simple_eval_block, _2D_simple_eval_block


log = logging.getLogger('westtools.w_kinavg')

from westtools.dtypes import iter_block_ci_dtype as ci_dtype

# From w_kinetics.
ed_list_dtype = numpy.dtype([('istate', numpy.uint16), ('fstate', numpy.uint16), ('duration', numpy.float64),
                             ('weight', numpy.float64), ('seg_id', seg_id_dtype)])
from westpa.binning import index_dtype
from westpa.kinetics._kinetics import _fast_transition_state_copy #@UnresolvedImport
from westpa.kinetics import find_macrostate_transitions

# From w_stateprobs
from westpa.binning import accumulate_state_populations_from_labeled


class DirectSubcommands(WESTSubcommand):
    '''Common argument processing for w_kinavg subcommands'''
    
    def __init__(self, parent):
        super(DirectSubcommands,self).__init__(parent)
        
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.progress = ProgressIndicatorComponent()
        
        self.output_filename = None
        # This is specific to the old w_kinavg, and isn't used by Kinetics.
        # This is actually applicable to both.
        self.assignment_filename = None
        
        self.output_file = None
        self.assignments_file = None
        
        self.evolution_mode = None
        
        self.mcbs_alpha = None
        self.mcbs_acalpha = None
        self.mcbs_nsets = None

        # Now we're adding in things that come from the old w_kinetics
        self.do_compression = True

        # We're going to try and re-use a lot of this code for the matrix commands, too, so.
        
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
        # We can do this with the output file, too...
        iogroup.add_argument('-k', '--kinetics', default=self.default_kinetics_file,
                            help='''Populations and transition rates are stored in KINETICS
                            (default: %(default)s).''')
        iogroup.add_argument('-o', '--output', dest='output', default=self.default_output_file,
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
        agroup = parser.add_argument_group('other options')
        agroup.add_argument('--config-from-file', dest='config_from_file', action='store_true', 
                            help='''Load bins/macrostates from a scheme specified in west.cfg.''')
        agroup.add_argument('--scheme-name', dest='scheme',
                            help='''Name of scheme specified in west.cfg.''')
        
    def open_files(self):
        #self.output_file = h5io.WESTPAH5File(self.output_filename, 'w', creating_program=True)
        self.output_file = h5io.WESTPAH5File(self.output_filename, 'a', creating_program=True)
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
        self.output_filename = os.path.join(path, 'kinavg.h5')
        self.kinetics_filename = os.path.join(path, 'kintrace.h5')
        self.assignments_filename = os.path.join(path, 'assign.h5')
        w_kinavg_config = { 'mcbs_alpha': 0.05, 'mcbs_nsets': 1000, 'evolution': 'cumulative', 'evol_window_frac': 1, 'step_iter': 1, 'bootstrap': True , 'correl': True, 'display_averages': False}
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
        self.display_averages = w_kinavg_config['display_averages']

def generate_future(work_manager, name, eval_block, kwargs):
    submit_kwargs = {'name': name}
    submit_kwargs.update(kwargs)
    future = work_manager.submit(eval_block, kwargs=submit_kwargs)
    return future

# Each of these blocks is responsible for submitting a set of calculations to be bootstrapped over for a particular type of calculation.
# A property which wishes to be calculated should adhere to this format.

def _rate_eval_block(iblock, start, stop, nstates, data_input, name, mcbs_alpha, mcbs_nsets, mcbs_acalpha, correl, **extra):
    # Our rate estimator is a little more complex, so we've defined a custom evaluation block for it,
    # instead of just using the block evalutors that we've imported.
    results = []
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            if istate == jstate: continue
            kwargs = { 'istate' : istate, 'jstate': jstate }
            dataset = {'dataset': data_input['dataset'][:, istate, jstate], 'pops': data_input['pops'] }
            ci_res = mcbs_ci_correl_rw(dataset,estimator=sequence_macro_flux_to_rate_bs,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=numpy.mean, pre_calculated=data_input['rates'][:,istate,jstate][numpy.isfinite(data_input['rates'][:,istate,jstate])], correl=correl, **kwargs)
            results.append((name, iblock, istate, jstate, (start,stop) + ci_res))

    return results


class DKinetics(DirectSubcommands):
    subcommand='kinetics'
    default_kinetics_file = 'kintrace.h5'
    default_output_file = 'kintrace.h5'
    help_text = 'calculate state-to-state kinetics by tracing trajectories'
    description = '''\
Calculate state-to-state rates and transition event durations by tracing
trajectories.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).

The output generated by this program is used as input for the ``w_kinavg``
tool, which converts the flux data in the output file into average rates
with confidence intervals. See ``w_kinavg trace --help`` for more 
information.

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "kintrace.h5") contains the
following datasets:

  ``/conditional_fluxes`` [iteration][state][state]
    *(Floating-point)* Macrostate-to-macrostate fluxes. These are **not**
    normalized by the population of the initial macrostate.

  ``/conditional_arrivals`` [iteration][stateA][stateB]
    *(Integer)* Number of trajectories arriving at state *stateB* in a given
    iteration, given that they departed from *stateA*.
    
  ``/total_fluxes`` [iteration][state]
    *(Floating-point)* Total flux into a given macrostate.
    
  ``/arrivals`` [iteration][state]
    *(Integer)* Number of trajectories arriving at a given state in a given
    iteration, regardless of where they originated.

  ``/duration_count`` [iteration]
    *(Integer)* The number of event durations recorded in each iteration.
    
  ``/durations`` [iteration][event duration]
    *(Structured -- see below)*  Event durations for transition events ending
    during a given iteration. These are stored as follows:
      
      istate
        *(Integer)* Initial state of transition event.
      fstate
        *(Integer)* Final state of transition event.
      duration
        *(Floating-point)* Duration of transition, in units of tau.
      weight
        *(Floating-point)* Weight of trajectory at end of transition, **not**
        normalized by initial state population.

Because state-to-state fluxes stored in this file are not normalized by
initial macrostate population, they cannot be used as rates without further
processing. The ``w_kinavg`` command is used to perform this normalization
while taking statistical fluctuation and correlation into account. See 
``w_kinavg trace --help`` for more information.  Target fluxes (total flux
into a given state) require no such normalization.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''
    
    def __init__(self, parent):
        super(DKinetics,self).__init__(parent)

    def go(self):
        pi = self.progress.indicator
        pi.new_operation('Initializing')
        with pi:
            self.data_reader.open('r')
            self.open_files()
            nstates = self.assignments_file.attrs['nstates']
            start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
            iter_count = stop_iter - start_iter
            durations_ds = self.output_file.create_dataset('durations', 
                                                           shape=(iter_count,0), maxshape=(iter_count,None),
                                                           dtype=ed_list_dtype,
                                                           chunks=(1,15360) if self.do_compression else None,
                                                           shuffle=self.do_compression,
                                                           compression=9 if self.do_compression else None)
            durations_count_ds = self.output_file.create_dataset('duration_count',
                                                                 shape=(iter_count,), dtype=numpy.int_, shuffle=True,compression=9)
            cond_fluxes_ds = self.output_file.create_dataset('conditional_fluxes',
                                                              shape=(iter_count,nstates,nstates), dtype=weight_dtype,
                                                              chunks=(h5io.calc_chunksize((iter_count,nstates,nstates),weight_dtype)
                                                                      if self.do_compression else None),
                                                              shuffle=self.do_compression,
                                                              compression=9 if self.do_compression else None)
            total_fluxes_ds = self.output_file.create_dataset('total_fluxes',
                                                              shape=(iter_count,nstates), dtype=weight_dtype,
                                                              chunks=(h5io.calc_chunksize((iter_count,nstates),weight_dtype)
                                                                      if self.do_compression else None),
                                                              shuffle=self.do_compression,
                                                              compression=9 if self.do_compression else None)

            cond_arrival_counts_ds = self.output_file.create_dataset('conditional_arrivals',
                                                                     shape=(iter_count,nstates,nstates), dtype=numpy.uint,
                                                                     chunks=(h5io.calc_chunksize((iter_count,nstates,nstates),
                                                                                                 numpy.uint)
                                                                             if self.do_compression else None),
                                                              shuffle=self.do_compression,
                                                              compression=9 if self.do_compression else None) 
            arrival_counts_ds = self.output_file.create_dataset('arrivals',
                                                                shape=(iter_count,nstates), dtype=numpy.uint,
                                                                chunks=(h5io.calc_chunksize((iter_count,nstates),
                                                                                            numpy.uint)
                                                                        if self.do_compression else None),
                                                                shuffle=self.do_compression,
                                                                compression=9 if self.do_compression else None)

            # copy state labels for convenience
            self.output_file['state_labels'] = self.assignments_file['state_labels'][...]

            # Put nice labels on things
            for ds in (self.output_file, durations_count_ds, cond_fluxes_ds, total_fluxes_ds):
                h5io.stamp_iter_range(ds, start_iter, stop_iter)

            # Calculate instantaneous rate matrices and trace trajectories
            last_state = None
            pi.new_operation('Tracing trajectories', iter_count)
            for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
                # Get data from the main HDF5 file
                iter_group = self.data_reader.get_iter_group(n_iter)
                seg_index = iter_group['seg_index']
                nsegs, npts = iter_group['pcoord'].shape[0:2] 
                weights = seg_index['weight']
                #parent_ids = seg_index['parent_id']
                parent_ids = self.data_reader.parent_id_dsspec.get_iter_data(n_iter)
                
                # Get bin and traj. ensemble assignments from the previously-generated assignments file
                assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
                bin_assignments = numpy.require(self.assignments_file['assignments'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                                dtype=index_dtype)
                label_assignments = numpy.require(self.assignments_file['trajlabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                                  dtype=index_dtype)
                state_assignments = numpy.require(self.assignments_file['statelabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                                  dtype=index_dtype)
                
                # Prepare to run analysis
                cond_fluxes = numpy.zeros((nstates,nstates), weight_dtype)
                total_fluxes = numpy.zeros((nstates,), weight_dtype)
                cond_counts = numpy.zeros((nstates,nstates), numpy.uint)
                total_counts = numpy.zeros((nstates,), numpy.uint)
                durations = []
    
                # Estimate macrostate fluxes and calculate event durations using trajectory tracing
                # state is opaque to the find_macrostate_transitions function            
                state = _fast_transition_state_copy(iiter, nstates, parent_ids, last_state)
                find_macrostate_transitions(nstates, weights, label_assignments, state_assignments, 1.0/(npts-1), state,
                                            cond_fluxes, cond_counts, total_fluxes, total_counts, durations)
                last_state = state
                
                # Store trace-based kinetics data
                cond_fluxes_ds[iiter] = cond_fluxes
                total_fluxes_ds[iiter] = total_fluxes
                arrival_counts_ds[iiter] = total_counts
                cond_arrival_counts_ds[iiter] = cond_counts
                
                durations_count_ds[iiter] = len(durations)
                if len(durations) > 0:
                    durations_ds.resize((iter_count, max(len(durations), durations_ds.shape[1])))
                    durations_ds[iiter,:len(durations)] = durations
                        
                # Do a little manual clean-up to prevent memory explosion
                del iter_group, weights, parent_ids, bin_assignments, label_assignments, state, cond_fluxes, total_fluxes
                pi.progress += 1

class AverageCommands(DirectSubcommands):
    default_output_file = 'direct.h5'

    def __init__(self, parent):
        # Ideally, this is stuff general to all the calculations we want to perform.
        super(AverageCommands,self).__init__(parent)
        self.kinetics_filename = None
        self.kinetics_file = None

    def open_assignments(self):
        # This seems to be stuff we're going to be using a lot, so.
        self.nstates = self.assignments_file.attrs['nstates']
        self.nbins = self.assignments_file.attrs['nbins']
        self.state_labels = self.assignments_file['state_labels'][...]
        assert self.nstates == len(self.state_labels)
        self.start_iter, self.stop_iter, self.step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

    def run_calculation(self, pi, nstates, start_iter, stop_iter, step_iter, dataset, eval_block, name, dim):
        #pi = self.progress.indicator
        
        start_pts = range(start_iter, stop_iter, step_iter)
        # Our evolution dataset!
        if dim == 2:
            evolution_dataset = numpy.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
        elif dim == 1:
            evolution_dataset = numpy.zeros((len(start_pts), nstates), dtype=ci_dtype)
        else:
            # Temp.
            print("What's wrong?")

        # This is appropriate for bootstrapped quantities, I think.
        all_items = numpy.arange(1,len(start_pts)+1)
        bootstrap_length = 0.5*(len(start_pts)*(len(start_pts)+1)) - numpy.delete(all_items, numpy.arange(1, len(start_pts)+1, step_iter))
        #with pi:
        if True:
            pi.new_operation('Calculating {}'.format(name), bootstrap_length[0])

            futures = []
            for iblock, start in enumerate(start_pts):
                stop = min(start+step_iter, stop_iter)
                if self.evolution_mode == 'cumulative':
                    windowsize = int(self.evol_window_frac * (stop - start_iter))
                    block_start = max(start_iter, stop - windowsize)
                else: # self.evolution_mode == 'blocked'
                    block_start = start

                # Create a basic set of kwargs for this iteration slice.
                future_kwargs = dict(iblock=iblock, start=block_start, stop=stop,
                                     nstates=nstates,
                                     mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                     mcbs_acalpha=self.mcbs_acalpha,
                                     correl=self.correl,name=name,
                                     data_input={})

                # Slice up the datasets for this iteration slice.
                # We're assuming they're all h5io iter blocked datasets; it's up to the calling routine
                # to ensure this is true.
                for key, value in dataset.iteritems():
                    future_kwargs['data_input'][key] = value.iter_slice(block_start,stop)

                # We create a future object with the appropriate name, and then append it to the work manager.
                futures.append(generate_future(self.work_manager, name, eval_block, future_kwargs))

            # Now, we wait to get the result back; we'll store it in the result, and return it.
            for future in self.work_manager.as_completed(futures):
                pi.progress += iblock / step_iter
                future_result = future.get_result(discard=True)
                # print(future_result)

                if dim == 2:
                    for result in future_result:
                        name,iblock,istate,jstate,ci_result = result
                        evolution_dataset[iblock, istate, jstate] = ci_result
                elif dim == 1:
                    for result in future_result:
                        name,iblock,istate,ci_result = result
                        evolution_dataset[iblock, istate] = ci_result


            return evolution_dataset


            
class DKinAvg(AverageCommands):
    subcommand = 'average'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'kintrace.h5'

    def go(self):
        pi = self.progress.indicator

        # We're initializing the various datasets...
        if True:
            self.open_files()
            self.open_assignments()
            # Obviously, this is for the conditional and total fluxes.  This is really all we need to sort for this.
            #pi.new_operation('Reading data')
            cond_fluxes = h5io.IterBlockedDataset(self.kinetics_file['conditional_fluxes'])
            cond_fluxes.cache_data()
            total_fluxes = h5io.IterBlockedDataset(self.kinetics_file['total_fluxes'])
            

            # This is necessary for both color and state populations...
            # ... but we also need this for the kinetics calculations.
            pops = h5io.IterBlockedDataset(self.assignments_file['labeled_populations'])
            pops.cache_data()
            pops.data = pops.data.sum(axis=2)

            # Here, we're pre-generating the information needed...
            rates = h5io.IterBlockedDataset.empty_like(cond_fluxes)
            rates.data = sequence_macro_flux_to_rate(cond_fluxes.data, pops.data[:self.nstates,:self.nbins])

            # As the dataset, just send in the rates and such for now.

        submit_kwargs = dict(pi=pi, nstates=self.nstates, start_iter=self.start_iter, stop_iter=self.stop_iter, 
                             step_iter=self.step_iter)

        with pi:
            submit_kwargs['dataset'] = {'dataset': cond_fluxes, 'pops': pops, 'rates': rates}
            rate_evol = self.run_calculation(eval_block=_rate_eval_block, name='Rate Evolution', dim=2, **submit_kwargs)
            self.output_file.replace_dataset('rate_evolution', data=rate_evol, shuffle=True, compression=9)

            submit_kwargs['dataset'] = {'dataset': cond_fluxes }
            rate_evol = self.run_calculation(eval_block=_2D_simple_eval_block, name='Conditional Flux Evolution', dim=2, **submit_kwargs)
            self.output_file.replace_dataset('conditional_flux_evolution', data=rate_evol, shuffle=True, compression=9)

            submit_kwargs['dataset'] = {'dataset': total_fluxes }
            rate_evol = self.run_calculation(eval_block=_1D_simple_eval_block, name='Target Flux Evolution', dim=1, **submit_kwargs)
            self.output_file.replace_dataset('target_flux_evolution', data=rate_evol, shuffle=True, compression=9)

class DStateProbs(AverageCommands):
    subcommand = 'stateprobs'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'kintrace.h5'

    def calculate_state_populations(self, pops):
            # ... but then this is how the state populations are done.
            # This was taken, more or less, from the old w_stateprobs
            iter_count = self.stop_iter-self.start_iter
            all_state_pops = numpy.empty((iter_count,self.nstates+1), weight_dtype)
            iter_state_pops = numpy.empty((self.nstates+1,), weight_dtype)
            avg_state_pops = numpy.zeros((self.nstates+1,), weight_dtype)
            pops.cache_data(max_size='available')
            state_map = self.assignments_file['state_map'][...]
            try:
                for iiter,n_iter in enumerate(xrange(self.start_iter,self.stop_iter)):
                    iter_state_pops.fill(0)
                    labeled_pops = pops.iter_entry(n_iter)
                    accumulate_state_populations_from_labeled(labeled_pops, state_map, iter_state_pops, check_state_map=False)
                    all_state_pops[iiter] = iter_state_pops
                    avg_state_pops += iter_state_pops
                    del labeled_pops
            finally:
                pops.drop_cache()

            state_pops = h5io.IterBlockedDataset.empty_like(pops)
            state_pops.data = all_state_pops
            return state_pops

    def go(self):
        pi = self.progress.indicator
        if True:
            self.open_files()
            self.open_assignments()
            # So far, we definitely need this boilerplate...
            #pi.new_operation('Reading data')

            # This is necessary for both color and state populations...
            pops = h5io.IterBlockedDataset(self.assignments_file['labeled_populations'])

            state_pops = self.calculate_state_populations(pops)

            # This now sorts it for the color populations
            pops.cache_data()
            pops.data = pops.data.sum(axis=2)

        submit_kwargs = dict(pi=pi,nstates=self.nstates, start_iter=self.start_iter, stop_iter=self.stop_iter, 
                             step_iter=self.step_iter, eval_block=_1D_simple_eval_block)

        with pi:
            submit_kwargs['dataset'] = {'dataset': pops}
            pop_evol = self.run_calculation(name='Color Probability Evolution', dim=1, **submit_kwargs)
            self.output_file.replace_dataset('color_prob_evolution', data=pop_evol, shuffle=True, compression=9)

            submit_kwargs['dataset'] = {'dataset': state_pops}
            pop_evol = self.run_calculation(name='State Probability Evolution', dim=1, **submit_kwargs)
            self.output_file.replace_dataset(name='state_pop_evolution', data=pop_evol, shuffle=True, compression=9)

class WDirect(WESTMasterCommand, WESTParallelTool):
    prog='w_direct'
    #subcommands = [AvgTraceSubcommand,AvgMatrixSubcommand]
    subcommands = [DKinetics, DKinAvg, DStateProbs]
    subparsers_title = 'direct kinetics analysis schemes'
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
    WDirect().main()
