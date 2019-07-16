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

from westtools import (WESTToolComponent, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)

import mclib
from mclib import mcbs_correltime, mcbs_ci_correl, _1D_simple_eval_block, _2D_simple_eval_block
from westpa import h5io

import numpy
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

# A function to just help with creating future objects for the work manager.

def generate_future(work_manager, name, eval_block, kwargs):
    submit_kwargs = {'name': name}
    submit_kwargs.update(kwargs)
    future = work_manager.submit(eval_block, kwargs=submit_kwargs)
    return future

    
class WESTKineticsBase(WESTSubcommand):
    '''
    Common argument processing for w_direct/w_reweight subcommands.
    Mostly limited to handling input and output from w_assign.
    '''
    
    def __init__(self, parent):
        super(WESTKineticsBase,self).__init__(parent)
        
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.progress = ProgressIndicatorComponent()
        
        self.output_filename = None
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
        
            
    def add_args(self, parser):
        self.progress.add_args(parser)
        self.data_reader.add_args(parser)
        self.iter_range.include_args['iter_step'] = True
        self.iter_range.add_args(parser)

        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')
        
        iogroup.add_argument('-o', '--output', dest='output', default=self.default_output_file,
                            help='''Store results in OUTPUT (default: %(default)s).''')

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

# This provides some convenience functions, modified from w_kinavg, to help with calculating evolution and averages for observables with the mclib library in a consistent manner.
# It's used in both w_direct and w_reweight.
class AverageCommands(WESTKineticsBase):
    default_output_file = 'direct.h5'

    def __init__(self, parent):
        # Ideally, this is stuff general to all the calculations we want to perform.
        super(AverageCommands,self).__init__(parent)
        self.kinetics_filename = None
        self.kinetics_file = None

    def add_args(self, parser):
        
        iogroup = parser.add_argument_group('input/output options')
        # self.default_kinetics_file will be picked up as a class attribute from the appropriate subclass        
        # We can do this with the output file, too...
        # ... by default, however, we're going to use {direct/reweight}.h5 for everything.
        # Modules which are called with different default values will, of course, still use those.
        iogroup.add_argument('-k', '--kinetics', default=self.default_kinetics_file,
                            help='''Populations and transition rates are stored in KINETICS
                            (default: %(default)s).''')

        cgroup = parser.add_argument_group('confidence interval calculation options')
        cgroup.add_argument('--disable-bootstrap', '-db', dest='bootstrap', action='store_const', const=False,
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
    
    def process_args(self, args):
        self.kinetics_filename = args.kinetics
                
        # Disable the bootstrap or the correlation analysis.
        self.mcbs_enable = args.bootstrap if args.bootstrap is not None else True
        self.do_correl = args.correl if args.correl is not None else True
        self.mcbs_alpha = args.alpha
        self.mcbs_acalpha = args.acalpha if args.acalpha else self.mcbs_alpha
        self.mcbs_nsets = args.nsets if args.nsets else mclib.get_bssize(self.mcbs_alpha)

        self.display_averages = args.display_averages
        
        self.evolution_mode = args.evolution_mode
        self.evol_window_frac = args.window_frac
        if self.evol_window_frac <= 0 or self.evol_window_frac > 1:
            raise ValueError('Parameter error -- fractional window defined by --window-frac must be in (0,1]')

    def stamp_mcbs_info(self, dataset):
        dataset.attrs['mcbs_alpha'] = self.mcbs_alpha
        dataset.attrs['mcbs_acalpha'] = self.mcbs_acalpha
        dataset.attrs['mcbs_nsets'] = self.mcbs_nsets
        
    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_filename, 'a', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.assignments_file = h5io.WESTPAH5File(self.assignments_filename, 'r')#, driver='core', backing_store=False)
        self.kinetics_file = h5io.WESTPAH5File(self.kinetics_filename, 'r')#, driver='core', backing_store=False)
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments data do not span the requested iterations')


    def open_assignments(self):
        # Actually, I should rename this, as we're not OPENING assignments.
        # This seems to be stuff we're going to be using a lot, so.
        self.nstates = self.assignments_file.attrs['nstates']
        self.nbins = self.assignments_file.attrs['nbins']
        self.state_labels = self.assignments_file['state_labels'][...]
        assert self.nstates == len(self.state_labels)
        self.start_iter, self.stop_iter, self.step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step

        # Import for the reweighting code.

        state_map = self.assignments_file['state_map'][...]
        # We've moved this into a different step so that it's compatible with
        # loading up from the all command.
        # Otherwise, we try to load the kinetics (since we're just mixing subclasses)
        # before it's actually run, and so we fail out.
        if not self.iter_range.check_data_iter_range_least(self.kinetics_file):
            raise ValueError('kinetics data do not span the requested iterations')

    def print_averages(self, dataset, header, dim=1):
        print(header)
        maxlabellen = max(list(map(len,self.state_labels)))
        for istate in range(self.nstates):
            if dim == 1:
                print('{:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e}) * tau^-1'
                        .format(self.state_labels[istate],
                        dataset['expected'][istate],
                        dataset['ci_lbound'][istate],
                        dataset['ci_ubound'][istate],
                        maxlabellen=maxlabellen))

            else:
                for jstate in range(self.nstates):
                    if istate == jstate: continue
                    print('{:{maxlabellen}s} -> {:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e}) * tau^-1'
                        .format(self.state_labels[istate], self.state_labels[jstate],
                        dataset['expected'][istate,jstate],
                        dataset['ci_lbound'][istate,jstate],
                        dataset['ci_ubound'][istate,jstate],
                        maxlabellen=maxlabellen))

    def run_calculation(self, pi, nstates, start_iter, stop_iter, step_iter, dataset, eval_block, name, dim, do_averages=False, **extra):
        
        # We want to use the same codepath to run a quick average as we do the longer evolution sets, so...
        if do_averages:
            start_pts = [start_iter, stop_iter]
        else:
            start_pts = list(range(start_iter, stop_iter, step_iter))
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
        if True:
            pi.new_operation('Calculating {}'.format(name), bootstrap_length[0])

            futures = []
            for iblock, start in enumerate(start_pts):
                stop = min(start+step_iter, stop_iter)
                if self.evolution_mode == 'cumulative' or do_averages == True:
                    windowsize = int(self.evol_window_frac * (stop - start_iter))
                    block_start = max(start_iter, stop - windowsize)
                else: # self.evolution_mode == 'blocked'
                    block_start = start

                # Create a basic set of kwargs for this iteration slice.
                future_kwargs = dict(iblock=iblock, start=block_start, stop=stop,
                                     nstates=nstates,
                                     mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                     mcbs_acalpha=self.mcbs_acalpha,
                                     do_correl=self.do_correl,name=name,
                                     mcbs_enable=self.mcbs_enable,
                                     data_input={},
                                     **extra)

                # Slice up the datasets for this iteration slice.
                # We're assuming they're all h5io iter blocked datasets; it's up to the calling routine
                # to ensure this is true.

                # Actually, I'm less sure how to handle this for pre-calculated datasets.  Need to consider this.  But for now...
                for key, value in dataset.items():
                    try:
                        future_kwargs['data_input'][key] = value.iter_slice(block_start,stop) if hasattr(value, 'iter_slice') else value[block_start:stop]
                    except:
                        future_kwargs['data_input'][key] = value.iter_slice(block_start,stop) if hasattr(value, 'iter_slice') else value[block_start:stop,:]
                    #print(future_kwargs['data_input'][key])

                # We create a future object with the appropriate name, and then append it to the work manager.
                futures.append(generate_future(self.work_manager, name, eval_block, future_kwargs))

            # Now, we wait to get the result back; we'll store it in the result, and return it.
            for future in self.work_manager.as_completed(futures):
                pi.progress += iblock / step_iter
                future_result = future.get_result(discard=True)

                if dim == 2:
                    for result in future_result:
                        name,iblock,istate,jstate,ci_result = result
                        evolution_dataset[iblock, istate, jstate] = ci_result
                elif dim == 1:
                    for result in future_result:
                        name,iblock,istate,ci_result = result
                        evolution_dataset[iblock, istate] = ci_result

        return evolution_dataset
