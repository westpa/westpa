from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTMasterCommand, WESTParallelTool, WESTDataReader, IterRangeSelection, WESTSubcommand
import sys, math
import numpy
from westpa import h5io

import westpa
from west.data_manager import weight_dtype, n_iter_dtype
import mclib
 
from westpa.kinetics import labeled_flux_to_rate, sequence_macro_flux_to_rate
from westpa.kinetics.matrates import get_macrostate_rates
from mclib import mcbs_correltime, mcbs_ci_correl


log = logging.getLogger('westtools.w_kinavg')

ci_dtype = numpy.dtype([('iter_start', n_iter_dtype),
                        ('iter_stop', n_iter_dtype),
                        ('expected', numpy.float64),
                        ('ci_lbound', numpy.float64),
                        ('ci_ubound', numpy.float64),
                        ('corr_len', n_iter_dtype)])

class KinAvgSubcommands(WESTSubcommand):
    '''Common argument processing for w_kinavg subcommands'''
    
    def __init__(self, parent):
        super(KinAvgSubcommands,self).__init__(parent)
        
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection() 
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
        cogroup.add_argument('--evolution-mode', choices=['cumulative', 'blocked', 'none'], default='none',
                             help='''How to calculate time evolution of rate estimates.
                             ``cumulative`` evaluates rates over windows starting with --start-iter and getting progressively
                             wider to --stop-iter by steps of --step-iter.
                             ``blocked`` evaluates rates over windows of width --step-iter, the first of which begins at
                             --start-iter.
                             ``none`` (the default) disables calculation of the time evolution of rate estimates.''')
        

    
    def process_args(self, args):
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        self.kinetics_file = h5io.WESTPAH5File(args.kinetics, 'r')
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args, default_iter_step=None)
        if self.iter_range.iter_step is None:
            #use about 10 blocks by default
            self.iter_range.iter_step = max(1, (self.iter_range.iter_stop - self.iter_range.iter_start) // 10)
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments data do not span the requested iterations')

        if not self.iter_range.check_data_iter_range_least(self.kinetics_file):
            raise ValueError('kinetics data do not span the requested iterations')
        
        self.mcbs_alpha = args.alpha
        self.mcbs_acalpha = args.acalpha if args.acalpha else self.mcbs_alpha
        self.mcbs_nsets = args.nsets if args.nsets else mclib.get_bssize(self.mcbs_alpha)
        
        self.evolution_mode = args.evolution_mode
        
def _eval_block(iblock, start, stop, nstates, rates, mcbs_alpha, mcbs_nsets, mcbs_acalpha):
    results = []
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            if istate == jstate: continue
            
            ci_res = mcbs_ci_correl(rates[:,istate,jstate],estimator=numpy.mean,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=numpy.mean)
            results.append((iblock, istate, jstate, (start,stop) + ci_res))
    return results
        
class AvgTraceSubcommand(KinAvgSubcommands):
    subcommand = 'trace'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'kintrace.h5'
    
    def __init__(self, parent):
        super(AvgTraceSubcommand,self).__init__(parent)
                        
    def go(self):
        nstates = self.assignments_file.attrs['nstates']
        state_labels = self.assignments_file['state_labels'][...]
        assert nstates == len(state_labels)
        start_iter, stop_iter, step_iter = self.iter_range.iter_start, self.iter_range.iter_stop, self.iter_range.iter_step
        
        fluxes = h5io.IterBlockedDataset(self.kinetics_file['trace_macro_fluxes'])
        fluxes.cache_data()
        pops = h5io.IterBlockedDataset(self.assignments_file['labeled_populations'])
        pops.cache_data()
        pops.data = pops.data.sum(axis=2)
        
        rates = h5io.IterBlockedDataset.empty_like(fluxes)
        rates.data = sequence_macro_flux_to_rate(fluxes.data, pops.data)
        
        avg_rates = numpy.zeros((nstates,nstates), dtype=ci_dtype)
        
        # Calculate overall average rates
        print('evaluating overall averages...')
        for istate in xrange(nstates):
            for jstate in xrange(nstates):
                if istate == jstate: continue
                
                ci_res = mcbs_ci_correl(rates.iter_slice(start_iter,stop_iter)[:,istate,jstate],estimator=numpy.mean,
                                        alpha=self.mcbs_alpha,n_sets=self.mcbs_nsets,autocorrel_alpha=self.mcbs_acalpha,
                                        subsample=numpy.mean)
                
                avg_rates[istate, jstate] = (start_iter, stop_iter) + ci_res
        self.output_file['avg_rates'] = avg_rates
        self.stamp_mcbs_info(self.output_file['avg_rates'])
        maxlabellen = max(map(len,state_labels))
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
        
        print('\nevaluating CI evolution...')
        
        start_pts = range(start_iter, stop_iter, step_iter)
        evol = numpy.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
        futures = []
        for iblock, start in enumerate(start_pts):
            if self.evolution_mode == 'cumulative':
                block_start = start_iter
            else: # self.evolution_mode == 'blocked'
                block_start = start
            stop = min(start+step_iter, stop_iter)
            
            future = self.work_manager.submit(_eval_block, kwargs=dict(iblock=iblock, start=block_start, stop=stop,
                                                                       nstates=nstates,
                                                                       rates=rates.iter_slice(block_start,stop),
                                                                       mcbs_alpha=self.mcbs_alpha, mcbs_nsets=self.mcbs_nsets,
                                                                       mcbs_acalpha=self.mcbs_acalpha))
            futures.append(future)
        
        for iresult, future in enumerate(self.work_manager.as_completed(futures)):
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\r{} of {} blocks done...'.format(iresult+1,len(start_pts)), end='')
            results = future.get_result(discard=True)
            for result in results:
                iblock, istate, jstate, ci_result = result 
                evol[iblock, istate, jstate] = ci_result
        if sys.stdout.isatty() and not westpa.rc.quiet_mode:
            print()
                    
        self.output_file.create_dataset('rate_evolution', data=evol, shuffle=True, compression=9)
        self.stamp_mcbs_info(self.output_file['rate_evolution'])

class WKinAvg(WESTMasterCommand, WESTParallelTool):
    prog='w_kinavg'
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
    
Each element in these datasets reflects a state-to-state rate matrix evaluted
over a certain window of iterations. The structure of the data is as follows:

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
    *alpha=0.5* corresponds to a 95% confidence interval.)

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
