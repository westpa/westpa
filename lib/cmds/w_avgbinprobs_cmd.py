from __future__ import division, print_function
import os, sys, argparse
from math import ceil, floor, log10
import numpy

import logging
log = logging.getLogger('w_avgbinprobs')

import wemd

parser = wemd.rc.common_arg_parser(description = '''\
Calculate the average probability distribution across bins and associated statistical error
''')
parser.add_argument('-o', '--output', dest='output_file', type=argparse.FileType('wt'), default=sys.stdout,
                    help='Store per-iteration output in OUTPUT_FILE (default: write to standard output).')
parser.add_argument('-l', '--labels', dest='print_labels', action='store_true',
                    help='Print bin labels (to stdout) corresponding to each column in the output '
                        +'(requires system driver to be available; default: do not print labels)')
parser.add_argument('-p', '--precision', dest='precision', type=int, 
                    help='Number of significant figures for probability display (default: 6)',
                    default=6)
parser.add_argument('-b', '--begin', '--start', dest='start_iter', type=int, default=1,
                    help='Begin averaging at iteration START_ITER (default: first iteration)')
parser.add_argument('-e', '--end', '--stop', dest='stop_iter', type=int,
                    help='Average through iteration STOP_ITER (default: last completed iteration)')
parser.add_argument('--confidence', dest='confidence', type=float, default=0.95,
                    help='Construct a confidence interval of width CONFIDENCE (default: 0.95=95%%)')
parser.add_argument('--bssize', dest='bssize', type=int,
                    help='Use a bootstrap of BSSIZE samples to calculate error (default: chosen from confidence)')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
args = parser.parse_args()

wemd.rc.config_logging(args, 'w_avgbinprobs')
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()

start_iter = args.start_iter
stop_iter = sim_manager.data_manager.current_iteration - 1
n_iters = stop_iter - start_iter + 1
if args.stop_iter and args.stop_iter <= stop_iter:
    stop_iter = args.stop_iter

sim_manager.load_system_driver()
bins = sim_manager.system.region_set.get_all_bins()
n_bins = len(bins)
        
confidence = args.confidence
alpha = 1 - confidence
if args.bssize:
    bssize = args.bssize
else:
    bssize = int(10**(ceil(-log10(alpha)) + 1))
sys.stdout.write('{:g}% confidence interval requested\n'.format(confidence*100))
sys.stdout.write('using bootstrap of {:d} samples to estimate confidence interval\n'.format(bssize))

# Retrieve bin probabilities
binprobs = numpy.empty((n_iters,n_bins), numpy.float64)
for n_iter in xrange(start_iter, stop_iter+1): 
    binprobs[n_iter-1] = sim_manager.data_manager.get_iter_group(n_iter)['bin_probs']

# Average over time
avg_binprobs = numpy.mean(binprobs, axis=0)
sumavg = avg_binprobs.sum()
assert len(avg_binprobs) == n_bins
sys.stdout.write('total average probability: {} (error: {})\n'.format(sumavg, abs(sumavg-1.0)))

# Use bootstrapping to obtain synthetic probability as a function of time, then average to obtain
# error bars on the average probability
syn_avgbinprobs = numpy.empty((bssize, n_bins), numpy.float64)
for ibin in xrange(0, n_bins):
    for isynth in xrange(0, bssize):
        # for each bin, select bssize samples of length n_iters, with replacement, from the original data
        indices = numpy.random.randint(n_iters, size=(n_iters,))
        syn_probs = binprobs[indices, ibin]
        syn_avgbinprobs[isynth,ibin] = syn_probs.mean() 
    syn_avgbinprobs[:,ibin].sort()

lbi = int(floor(bssize*alpha/2))
ubi = int(ceil(bssize*(1-alpha/2)))

lb_binprobs = syn_avgbinprobs[lbi,:]
ub_binprobs = syn_avgbinprobs[ubi,:]
    
if not args.suppress_headers:
    args.output_file.write('''\
# WE bin population analysis
# Iterations {start_iter} -- {stop_iter} (inclusive)
# Confidence level: {confidence}
# Number of bootstrap data sets: {bssize}
# ----
'''.format(start_iter = start_iter, stop_iter = stop_iter, confidence = confidence, bssize = bssize))

    if args.print_labels:
        maxwidth = int(ceil(log10(n_bins)))
        sys.stdout.write('# bin labels:\n')
        for (ibin, bin) in enumerate(bins):
            args.output_file.write('#  bin {0:<{maxwidth}d}: {1!s}\n'.format(ibin, bin.label, maxwidth=maxwidth))
        args.output_file.write('# ----\n')
    args.output_file.write('''\
# column 0:  bin index
# column 1:  average population
# column 2:  lower bound of confidence interval
# column 3:  upper bound of confidence interval
# column 4:  width of confidence interval
# column 5:  relative width of confidence interval (width/average)
# column 6:  symmetrized error [max(upper - average, average-lower)]
''')
    
max_bin_width = int(ceil(log10(n_bins)))
prob_width = min(12, args.precision+5)
prob_fmt = '{:d}.{:d}e'.format(prob_width, args.precision)
for ibin in xrange(0, n_bins):
    avg = avg_binprobs[ibin]
    lb = lb_binprobs[ibin]
    ub = ub_binprobs[ibin]
    rel_width = (ub-lb)/avg if avg != 0 else 0.0
    args.output_file.write(( '{ibin:<{max_bin_width}d}    {avg:{prob_fmt}}    {lb:{prob_fmt}}    {ub:{prob_fmt}}    '
                            +'{ci_width:{prob_fmt}}    {rel_ci_width:{prob_fmt}}    {symm_error:{prob_fmt}}\n')
                            .format(ibin=ibin, avg=avg, lb=lb, ub=ub, ci_width=ub-lb, rel_ci_width=rel_width,
                                    symm_error=max(ub-avg, avg-lb),
                                    max_bin_width=max_bin_width, prob_fmt=prob_fmt))

