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

import numpy

from west.data_manager import weight_dtype
from westtools import (WESTParallelTool, WESTDataReader, IterRangeSelection, ProgressIndicatorComponent)
from westpa import h5io
from westpa.binning import accumulate_state_populations_from_labeled

import mclib
from mclib import mcbs_ci_correl_rw


log = logging.getLogger('westtools.w_stateprobs')

from westtools.dtypes import iter_block_ci_dtype as ci_dtype

class StateProbTool(WESTParallelTool):
    prog='w_stateprobs'
    description = '''\
Calculate average populations and associated errors in state populations from
weighted ensemble data. Bin assignments, including macrostate definitions,
are required. (See "w_assign --help" for more information).

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, usually "stateprobs.h5") contains the following
dataset:

  /avg_state_pops [state]
    (Structured -- see below) Population of each state across entire
    range specified.

If --evolution-mode is specified, then the following additional dataset is
available:

  /state_pop_evolution [window][state]
    (Structured -- see below). State populations based on windows of
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
    
    def __init__(self):
        super(StateProbTool,self).__init__()
        
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.progress = ProgressIndicatorComponent()
        
        self.output_filename = None
        self.kinetics_filename = None
        
        self.output_file = None
        self.assignments_file = None
        
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
        iogroup.add_argument('-o', '--output', dest='output', default='stateprobs.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

        
        cgroup = parser.add_argument_group('confidence interval calculation options')
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
        
    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_filename, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.assignments_file = h5io.WESTPAH5File(self.assignments_filename, 'r')#, driver='core', backing_store=False)
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments data do not span the requested iterations')

    
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

        self.mcbs_alpha = args.alpha
        self.mcbs_acalpha = args.acalpha if args.acalpha else self.mcbs_alpha
        self.mcbs_nsets = args.nsets if args.nsets else mclib.get_bssize(self.mcbs_alpha)
        
        self.evolution_mode = args.evolution_mode
        
    def calc_state_pops(self):
        start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop
        nstates = self.nstates
        state_map = self.state_map
        iter_count = stop_iter-start_iter
        
        pi = self.progress.indicator
        pi.new_operation('Calculating state populations')
        pops = h5io.IterBlockedDataset(self.assignments_file['labeled_populations'])
        
        iter_state_pops = numpy.empty((nstates+1,), weight_dtype)
        all_state_pops = numpy.empty((iter_count,nstates+1), weight_dtype)
        avg_state_pops = numpy.zeros((nstates+1,), weight_dtype)
        pops.cache_data(max_size='available')
        try:
            for iiter,n_iter in enumerate(xrange(start_iter,stop_iter)):
                iter_state_pops.fill(0)
                labeled_pops = pops.iter_entry(n_iter)
                accumulate_state_populations_from_labeled(labeled_pops, state_map, iter_state_pops, check_state_map=False)
                all_state_pops[iiter] = iter_state_pops
                avg_state_pops += iter_state_pops
                del labeled_pops
                pi.progress += 1
        finally:
            pops.drop_cache()
        self.output_file.create_dataset('state_pops', data=all_state_pops, compression=9, shuffle=True)
        h5io.stamp_iter_range(self.output_file['state_pops'], start_iter, stop_iter)
        
        self.all_state_pops = all_state_pops
        avg_state_pops = numpy.zeros((nstates+1,), ci_dtype)
        pi.new_operation('Calculating overall average populations and CIs', nstates)
#        futures = []
#         for istate in xrange(nstates):
#             futures.append(self.work_manager.submit(_eval_block,kwargs=dict(iblock=None,istate=istate,
#                                                                             start=start_iter,stop=stop_iter,
#                                                                             state_pops=all_state_pops[:,istate],
#                                                                             mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
#                                                                             mcbs_acalpha = self.mcbs_acalpha)))
#         for future in self.work_manager.as_completed(futures):
        def taskgen():
            for istate in xrange(nstates):
                yield (_eval_block, (), dict(iblock=None,istate=istate,
                                             start=start_iter,stop=stop_iter,
                                             state_pops=all_state_pops[:,istate],
                                             mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                             mcbs_acalpha = self.mcbs_acalpha))
        for future in self.work_manager.submit_as_completed(taskgen(), self.max_queue_len):
            (_iblock,istate,ci_res) = future.get_result(discard=True)
            avg_state_pops[istate] = ci_res
            pi.progress += 1
        self.output_file['avg_state_pops'] = avg_state_pops
        self.stamp_mcbs_info(self.output_file['avg_state_pops'])
        pi.clear()
        
        maxlabellen = max(map(len,self.state_labels))
        print('average state populations:')
        for istate in xrange(nstates):
            print('{:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e})'
                  .format(self.state_labels[istate],
                          avg_state_pops['expected'][istate],
                          avg_state_pops['ci_lbound'][istate],
                          avg_state_pops['ci_ubound'][istate],
                          maxlabellen=maxlabellen))
        
    def calc_evolution(self):
        nstates = self.nstates
        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step
        start_pts = range(start_iter, stop_iter, step_iter)

        pop_evol = numpy.zeros((len(start_pts), nstates), dtype=ci_dtype)

        pi = self.progress.indicator
        pi.new_operation('Calculating population evolution', len(start_pts)*nstates)
#         futures = []
#         for iblock, start in enumerate(start_pts):
#             if self.evolution_mode == 'cumulative':
#                 block_start = start_iter
#             else: # self.evolution_mode == 'blocked'
#                 block_start = start
#             stop = min(start+step_iter, stop_iter)
# 
#             for istate in xrange(nstates):
#                 future = self.work_manager.submit(_eval_block,kwargs=dict(iblock=iblock,istate=istate,
#                                                                           start=block_start,stop=stop,
#                                                                           state_pops=self.all_state_pops[block_start-start_iter:stop-start_iter,istate],
#                                                                           mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
#                                                                           mcbs_acalpha = self.mcbs_acalpha))
#                 futures.append(future)
        def taskgen():
            for iblock, start in enumerate(start_pts):
                if self.evolution_mode == 'cumulative':
                    block_start = start_iter
                else: # self.evolution_mode == 'blocked'
                    block_start = start
                stop = min(start+step_iter, stop_iter)
     
                for istate in xrange(nstates):
                    yield (_eval_block,(),dict(iblock=iblock,istate=istate,
                                               start=block_start,stop=stop,
                                               state_pops=self.all_state_pops[block_start-start_iter:stop-start_iter,istate],
                                               mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                               mcbs_acalpha = self.mcbs_acalpha))
        #for future in self.work_manager.as_completed(futures):
        for future in self.work_manager.submit_as_completed(taskgen(), self.max_queue_len):
            (iblock,istate,ci_res) = future.get_result(discard=True)
            pop_evol[iblock,istate] =  ci_res
            pi.progress += 1

        self.output_file.create_dataset('state_pop_evolution', data=pop_evol, shuffle=True, compression=9)
        pi.clear()

    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            nstates = self.nstates = self.assignments_file.attrs['nstates']

            state_labels = self.state_labels = self.assignments_file['state_labels'][...]
            state_map = self.state_map = self.assignments_file['state_map'][...]
            if (state_map > nstates).any():
                raise ValueError('invalid state mapping')

            # copy metadata to output
            self.output_file.attrs['nstates'] = nstates
            self.output_file['state_labels'] = state_labels

            # calculate overall averages
            self.calc_state_pops()

            # calculate evolution, if requested
            if self.evolution_mode != 'none' and self.iter_range.iter_step:
                self.calc_evolution()

def _eval_block(iblock, istate, start, stop, state_pops, mcbs_alpha, mcbs_nsets, mcbs_acalpha):
    #ci_res = mcbs_ci_correl(state_pops,estimator=numpy.mean,alpha=mcbs_alpha,n_sets=mcbs_nsets,
    #                        autocorrel_alpha=mcbs_acalpha,subsample=numpy.mean)

    ci_res = mcbs_ci_correl({'dataset': state_pops},estimator=(lambda stride, dataset: numpy.mean(dataset)), ,alpha=mcbs_alpha,n_sets=mcbs_nsets,
                            autocorrel_alpha=mcbs_acalpha,subsample=numpy.mean, pre_calculated=state_pops, correl=False)
    return (iblock,istate,(start,stop)+ci_res)

if __name__ == '__main__':
    StateProbTool().main()
