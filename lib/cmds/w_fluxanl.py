'''
w_fluxanl: Obtain rate constants from flux into target states.

Monte Carlo bootstrapping is used to obtain confidence intervals on the
resulting rates.
'''

import os, sys, operator, argparse
import numpy
import wemd, wemdtools
from wemdtools.stats.mcbs import bootstrap_ci
from math import floor, ceil, log10
from itertools import izip

import logging
log = logging.getLogger('w_fluxanl')

parser = wemd.rc.common_arg_parser('w_fluxanl', description='''Obtain rate constants from flux into target states.''')
parser.add_argument('-o', '--output', dest='output_pattern', default='avgflux_%s.txt',
                    help='Store output in OUTPUT_PATTERN, which must contain a printf-style pattern '
                        +'which will be replaced by the name/index of the target state to which results correspond '
                        +' (default: %(default)s).')
parser.add_argument('--output-fluxes', dest='flux_output_pattern',
                    help='Store target flux and number of target hits for each iteration in FLUX_OUTPUT_PATTERN, '
                        +'which must contain a printf-style pattern which will be replaced by the name/index of a '
                        +'target state (default: no output).')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')

parser.add_argument('-t', '--tau', dest='tau', type=float, default=1.0,
                    help='Propagation/resampling timestep (default: 1.0).')
parser.add_argument('-b', '--begin', '--start', dest='start', type=int, default=1,
                    help='Begin flux averaging at iteration START (default: 1)')
parser.add_argument('-e', '--end', '--stop', dest='stop', type=int,
                    help='Stop flux averaging at iteration STOP (default: last completed iteration).')
parser.add_argument('-s', '--step', dest='step', type=int, default=1,
                    help='Report rate every STEP resampling steps (default: 1).')

wemdtools.stats.mcbs.add_mcbs_options(parser)

parser.add_argument('datafile', nargs='*',
                    help='Read flux information from DATAFILE(s) (default: load WEMD HDF5 file specified in wemd.cfg).')

args = parser.parse_args()

if args.datafile:
    from wemd.util.config_dict import ConfigDict
    runtime_config = ConfigDict()
    runtime_config['data.h5file'] = args.datafile[0]
    print('reading WEMD data from', args.datafile[0])
else:
    print('reading data from WEMD simulation')
    runtime_config = wemd.rc.read_config(args.run_config_file)
    
runtime_config.update_from_object(args)
wemd.rc.config_logging(args, 'w_fluxanl')
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()
data_manager = wemdtools.data_manager.CachingDataReader(sim_manager.data_manager)

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
max_iter_width = max(8,len(str(max_iter)))
    
n_iters = stop_iter - start_iter + 1
if n_iters == 0:
    sys.stdout.write('No data to analyze; exiting.  If this seems wrong, check --start and --stop.\n')
    sys.exit(0)
    
n_targets = len(data_manager.get_iter_group(start_iter)['recycling'])
fluxes = numpy.empty((n_iters,), numpy.float64)
counts = numpy.empty((n_iters,), numpy.uint)
for itarget in xrange(n_targets):
    for (iiter, n_iter) in enumerate(xrange(start_iter, stop_iter+1)):
        recycling = data_manager.get_iter_group(n_iter)['recycling'][itarget]
        fluxes[iiter] = recycling['weight']
        counts[iiter] = recycling['count']

    fluxes /= args.tau
    
    # Write per-iteration fluxes, if requested        
    if args.flux_output_pattern:
        flux_output = file(args.flux_output_pattern % itarget, 'wt')
        if not args.suppress_headers:
            flux_output.write('''\
# Weighted ensemble flux values
# Target index: {itarget:d}
# tau:          {tau:g}
# ----
# column 0: iteration
# column 1: flux
# column 2: count
# ----            
'''.format(itarget=itarget,tau=args.tau))
        for (iiter, (flux, count)) in enumerate(izip(fluxes,counts)):
            n_iter = iiter + start_iter
            flux_output.write('{n_iter:<{max_iter_width}d}    {flux:<24.16g}    {count:d}\n'.format(max_iter_width = max_iter_width,
                                                                                                    n_iter = long(n_iter),
                                                                                                    flux = float(flux), 
                                                                                                    count = long(count)))

        flux_output.close()

    # Perform Monte Carlo bootstrapping on fluxes if requested
    if args.output_pattern:
        stats_output = file(args.output_pattern % itarget, 'wt')
        
        if not args.suppress_headers:
            stats_output.write('''\
# Weighted ensemble flux values
# Target index:      {itarget:d}
# tau:               {tau:g}
# Initial iteration: {start_iter:d}
# Confidence level:  {confidence:g}%
# ----
# column 0: iteration
# column 1: average flux
# column 2: lower bound of CI
# column 3: upper bound of CI
# column 4: width of CI
# column 5: relative width of CI [abs(width/average)]
# column 6: symmetrized error [max(upper bound - average, average - lower bound)]
# ----            
'''.format(itarget=itarget, tau=args.tau, start_iter=start_iter,confidence=args.confidence*100))
    
    
        linefmt = '    '.join(['{n_iter:<8d}'] + ['{:24.16g}']*6 + ['\n'])
        for imaxiter in xrange(args.step, n_iters+args.step-1, args.step):
            imaxiter = min(imaxiter, n_iters - 1) 
            
            cidata = bootstrap_ci(numpy.mean, fluxes[:imaxiter], alpha, bssize, extended_output = True)
            
            stats_output.write(linefmt.format(*map(float,cidata), n_iter = start_iter + imaxiter))
        
        stats_output.close()

    

