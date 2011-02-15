'''
w_fluxanl: Obtain rate constants from flux into target states.

Monte Carlo bootstrapping is used to obtain confidence intervals on the
resulting rates.
'''

import os, sys, operator, argparse
import numpy
import wemd
from math import floor, ceil, log10
from itertools import izip

import logging
log = logging.getLogger('w_fluxanl')

from collections import namedtuple

parser = wemd.rc.common_arg_parser('w_fluxanl', description='''Obtain rate constants from flux into target states.''')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
parser.add_argument('--output-fluxes', dest='flux_output_file', type=argparse.FileType('wt'),
                    help='Store target flux and number of target hits for each iteration in FLUX_OUTPUT_FILE '
                        +'(default: no output).')
parser.add_argument('-t', '--tau', dest='tau', type=float, default=1.0,
                    help='Propagation/resampling timestep (default: 1.0).')
parser.add_argument('-b', '--begin', '--start', dest='start', type=int, default=1,
                    help='Begin flux averaging at iteration START (default: 1)')
parser.add_argument('-e', '--end', '--stop', dest='stop', type=int,
                    help='Stop flux averaging at iteration STOP (default: last completed iteration).')
parser.add_argument('-s', '--step', dest='step', type=int,
                    help='Report rate every STEP resampling steps (default: report only final estimate).')
parser.add_argument('--confidence', dest='confidence', type=float, default=0.95,
                    help='Construct a confidence interval of width CONFIDENCE (default: 0.95=95%%)')
parser.add_argument('--bssize', dest='bssize', type=int,
                    help='Use a bootstrap of BSSIZE samples to calculate error (default: chosen from confidence)')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
args = parser.parse_args()
output_file = args.output_file

wemd.rc.config_logging(args)
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()
h5file = sim_manager.data_manager.h5file

confidence = args.confidence
alpha = 1 - confidence
if args.bssize:
    bssize = args.bssize
else:
    bssize = int(10**(ceil(-log10(alpha)) + 1))
sys.stdout.write('tau: {:g}\n'.format(args.tau))
sys.stdout.write('{:g}% confidence interval requested\n'.format(confidence*100))
sys.stdout.write('using bootstrap of {:d} samples to estimate confidence interval\n'.format(bssize))

start_iter = args.start
max_iter = sim_manager.data_manager.current_iteration - 1
if args.stop:
    stop_iter = min(max_iter, args.stop)
else:
    stop_iter = max_iter

summary_table = h5file['/summary']
all_fluxes = summary_table['target_flux'][:max_iter] / args.tau
all_counts = summary_table['target_hits'][:max_iter]

max_iter_width = max(10, int(ceil(log10(max_iter))))
max_count_width = max(6, int(ceil(log10(all_counts.max()))))

if args.flux_output_file is not None:
    if not args.suppress_headers:
        args.flux_output_file.write('{:<{max_iter_width}s}  {:<24s}  {}\n'.format('#Iteration', 'Flux', 'Counts',
                                                                                max_iter_width = max_iter_width))
    for (iiter, (flux, count)) in enumerate(izip(all_fluxes, all_counts)):
        args.flux_output_file.write('{n_iter:<{max_iter_width}d}  {flux:<24.16g}  {count:>{max_count_width}d}\n'\
                                    .format(n_iter = iiter+1, flux = flux, count=long(count),
                                            max_iter_width=max_iter_width, max_count_width=max_count_width))

if args.step:
    ub_taus = range(start_iter, stop_iter, args.step)
    if ub_taus[-1] < stop_iter:
        ub_taus.append(stop_iter)
else:
    ub_taus = [stop_iter]

if not args.suppress_headers:
    args.output_file.write('''\
# WE flux analysis
# tau: {tau:g}
# first flux value considered: iteration {start_iter}
# confidence level: {confidence}
# number of bootstrap data sets: {bssize} 
# ----
# column 0:  last iteration considered
# column 1:  average flux
# column 2:  lower bound of confidence interval
# column 3:  upper bound of confidence interval
# column 4:  width of confidence interval
# column 5:  relative width of confidence interval (width/average)
# column 5:  symmetrized error [max(upper - average, average - lower)]
'''.format(tau=args.tau, start_iter=start_iter, confidence=confidence, bssize=bssize))
    
for ub_tau in ub_taus:
    fluxes = all_fluxes[start_iter-1:ub_tau]
    counts = all_counts[start_iter-1:ub_tau]
    
    mean_flux = fluxes.mean()
    bsmean = numpy.empty((bssize,), numpy.float64)
    for i in xrange(0, bssize):
        indices = numpy.random.randint(len(fluxes), size=(len(fluxes),))
        bsmean[i] = fluxes[indices].mean() 
    bsmean.sort()
    lb = bsmean[int(floor(bssize*alpha/2))]
    ub = bsmean[int(ceil(bssize*(1-alpha/2)))] 
    
    args.output_file.write(('{ub_tau:<{max_iter_width}d}  {avg_flux:10.6e}  {lb:10.6e}  {ub:10.6e}  '
                            +'{ci_width:10.6e}  {rel_ci_width:<10.6e}  {symm_error:<10.6e}\n')\
                           .format(start_iter = start_iter, ub_tau = ub_tau, avg_flux = mean_flux, lb = lb, ub = ub,
                                   ci_width = ub-lb, rel_ci_width = (ub-lb)/mean_flux, 
                                   symm_error=max(ub - mean_flux, mean_flux - lb),
                                   max_iter_width=max_iter_width))
    
