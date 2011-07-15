from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math
import numpy, h5py
import wemd, wemdtools
from wemdtools.transitions.transacc import TransitionEventAccumulator
from wemdtools.trajectories.trajtree import TrajTree
from wemdtools.stats.mcbs import bootstrap_ci 
from math import ceil, log10

import logging
log = logging.getLogger('w_ttimes')

def print_2d_int_array(array, width=6, sep='  ', file=sys.stdout):
    for row in array:
        file.write('{}\n'.format(sep.join('{:{width}d}'.format(long(field), width=width) for field in row)))
    try:
        file.flush()
    except AttributeError:
        pass

def print_2d_float_array(array, precision=6, sep='  ', file=sys.stdout):
    width = precision+6
    for row in array:
        file.write('{}\n'.format(sep.join('{:{width}.{precision}e}'.format(float(field),
                                                                           width = width,
                                                                           precision=precision) for field in row)))
    try:
        file.flush()
    except AttributeError:
        pass
            
def accumulate_transitions_segment(segment, children, history, transacc, data_manager, trajtree, binprobs, start_iter, args):
    n_visited = trajtree.segments_visited
    n_to_visit = trajtree.segments_to_visit
    if not args.quiet_mode and n_visited % 1000 == 0:
        nsegwidth = int(ceil(log10(n_to_visit)))
        pct_visited = n_visited / n_to_visit * 100
        if args.whole_only:
            sys.stdout.write('{:>{nsegwidth}d} of a maximum of {:<{nsegwidth}d} segments analyzed\n'
                                   .format(n_visited, n_to_visit, nsegwidth=nsegwidth))
        else:
            sys.stdout.write('{:>{nsegwidth}d} of {:<{nsegwidth}d} segments ({:5.1f}%) analyzed\n'
                                   .format(n_visited, n_to_visit, pct_visited, nsegwidth=nsegwidth))
        sys.stdout.flush()

    if len(history) == 0:
        transacc.accumulate_transitions(segment.pcoord, weight=segment.weight, 
                                        region_weights=binprobs[segment.n_iter-start_iter],
                                        continuation = False, n_iter = segment.n_iter)
    else:
        transacc.accumulate_transitions(segment.pcoord[1:], weight=segment.weight, 
                                        region_weights=binprobs[segment.n_iter-start_iter], 
                                        continuation = True, n_iter = segment.n_iter)
        
    # Delete the progress coordinate from processed segments so that memory use doesn't explode
    del segment.pcoord
    segment.pcoord = None


parser = wemd.rc.common_arg_parser('w_ttimes', description='''
Perform lifetime, transition, and kinetic analysis on WEMD data.''')
parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                    help='Do not emit status messages (default: emit status every 1000 WE segments)')

parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: analysis.h5)',
                    default='analysis.h5')
parser.add_argument('-g', '--output-group', dest='output_group',
                    help='Store output in the given HDF5 group (default: w_ttimes)',
                    default='w_ttimes')
parser.add_argument('-L', '--ltstats', dest='lifetime_stats', type=argparse.FileType('wt'), default='lifetimes.txt',
                    help='Store per-bin lifetime statistics in LIFETIME_STATS (default: lifetimes.txt).')
parser.add_argument('-E', '--edstats', dest='ed_stats', type=argparse.FileType('wt'), default='eds.txt',
                    help='Store per-bin event duration statistics in ED_STATS (default: eds.txt)')
parser.add_argument('-F', '--fluxstats', dest='flux_stats', type=argparse.FileType('wt'), default='fluxes.txt',
                    help='Store per-bin flux statistics in FLUX_STATS (default: fluxes.txt)')
parser.add_argument('-R', '--ratestats', dest='rate_stats', type=argparse.FileType('wt'), default='rates.txt',
                    help='Store per-bin rate statistics in RATE_STATS (default: rates.txt)')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to (text) output files (default: write headers).')
parser.add_argument('-l', '--labels', dest='print_labels', action='store_true',
                    help='Print bin labels in headers of (text) output files (default: do not print labels)')

wemdtools.bins.add_region_set_options(parser)

parser.add_argument('--start', dest='start', type=int, default=0,
                    help='Start at WEMD iteration START (default: the beginning)')
parser.add_argument('--stop', dest='stop', type=int,
                    help='Stop after WEMD iteration STOP (default: the last completed iteration)')
parser.add_argument('--whole-only', dest='whole_only', action='store_true',
                    help='Consider entire trajectories only (default: consider any trajectory '
                        +'alive at START)')

parser.add_argument('--cache-pcoords', dest='cache_pcoords', action='store_true',
                    help = '(WE only) cache progress coordinate data in memory for improved analysis speed; '
                         + 'generally requires as much free RAM as the size of the HDF5 file (default: do not cache')

parser.add_argument('--dt', dest='dt', type=float, default=1.0,
                    help = 'When displaying times, consider the sampling timestep (time between pcoord records) to be DT'
                         + ' (default: 1.0)')

wemdtools.stats.mcbs.add_mcbs_options(parser)

parser.add_argument('datafile', nargs='*',
                    help='Read progress coordinate from DATAFILE(s) (default: load WEMD HDF5 file specified in wemd.cfg).')
args = parser.parse_args()

wemd.rc.config_logging(args, 'w_ttimes')

runtime_config = None
sim_manager = None
region_set = wemdtools.bins.get_region_set_from_args(args, status_stream=sys.stderr)
n_bins = len(region_set.get_all_bins())

if args.datafile:
    from wemd.util.config_dict import ConfigDict
    if runtime_config is None: runtime_config = ConfigDict()
    runtime_config['data.h5file'] = args.datafile[0]
    print('reading WEMD data from', args.datafile[0])
else:
    print('reading data from WEMD simulation')
    if runtime_config is None:
        runtime_config = wemd.rc.read_config(args.run_config_file)
                    
if sim_manager is None:
    sim_manager = wemd.rc.load_sim_manager(runtime_config)

sim_manager.load_data_manager()
sim_manager.data_manager.open_backing(mode='r')

start_iter = args.start or 1
stop_iter = sim_manager.data_manager.current_iteration - 1 if not args.stop else args.stop    
stop_iter  = min(stop_iter, sim_manager.data_manager.current_iteration-1)

if start_iter == stop_iter:
    sys.stdout.write('start and stop iterations are the same -- no data to analyze')
    sys.exit(0)

if args.whole_only:
    sys.stdout.write('considering whole trajectories only\n')

if args.cache_pcoords:
    sys.stdout.write('will cache pcoord data in memory\n')
    

output_h5 = h5py.File(args.output_file)
if args.output_group in output_h5:
    sys.stdout.write('overwriting data in group {}\n'.format(args.output_group))
    del output_h5[args.output_group]
output_group = output_h5.create_group(args.output_group)

tree = TrajTree(sim_manager.data_manager, cache_pcoords = args.cache_pcoords)
transacc = TransitionEventAccumulator(region_set, output_group)

# Use the caching wrapper, which is conveniently created around the true data manager when instantiating TrajTree
data_manager = tree.data_manager

# Determine the shape of pcoords - (n_segs, pcoord_len, pcoord_ndim)
pcoord_ds = data_manager.get_pcoord_dataset(1)
pcoord_len = pcoord_ds.shape[1]
pcoord_ndim = pcoord_ds.shape[2]

if not args.quiet_mode:
    sys.stdout.write('determining bin probabilities\n')
binprobs = numpy.empty((stop_iter-start_iter+1, pcoord_len, n_bins),dtype=numpy.float64)

for n_iter in xrange(start_iter, stop_iter+1):
    if not args.quiet_mode and n_iter % 10 == 0:
        sys.stdout.write('iteration {}\n'.format(n_iter))
    i_iter = n_iter - start_iter
    pcoords = data_manager.get_pcoord_array(n_iter)
    weights = data_manager.get_seg_index(n_iter)['weight']
    for ti in xrange(0,pcoord_len):
        index_map = region_set.map_to_all_indices(pcoords[:,ti,:])
        for si in xrange(0, len(index_map)):
            binprobs[i_iter,ti,index_map[si]] += weights[si]
    
sys.stdout.write('tracing trajectories\n')
tree.trace_trajectories(start_iter, stop_iter, 
                        callable=accumulate_transitions_segment,
                        get_state=transacc.get_state,
                        set_state=transacc.set_state,
                        args=(transacc,data_manager,tree,binprobs,start_iter,args),
                        whole_only = bool(args.whole_only))
transacc.flush_transition_data()
data_manager.clear_cache()
sys.stdout.write('{} segments analyzed\n'.format(tree.segments_visited))

# Prepare output files
for (output_file,title, ndim) in ((args.lifetime_stats, 'lifetime', 1),
                            (args.ed_stats, 'event duration', 2),
                            (args.flux_stats, 'flux', 2),
                            (args.rate_stats, 'rate', 2)):
    if not args.suppress_headers:
        if ndim == 1:
            output_file.write('''\
# column 0:  bin index
# column 1:  mean {title}
# column 2:  lower bound of confidence interval
# column 3:  upper bound of confidence interval
# column 4:  width of confidence interval
# column 5:  relative width of confidence interval (width/average)
# column 6:  symmetrized error [2*max(upper - average, average-lower)]
'''.format(title=title))
        else:
            output_file.write('''\
# column 0:  initial bin index
# column 1:  final bin index
# column 2:  mean {title}
# column 3:  lower bound of confidence interval
# column 4:  upper bound of confidence interval
# column 5:  width of confidence interval
# column 6:  relative width of confidence interval (width/average)
# column 7:  symmetrized error [2*max(upper - average, average-lower)]
'''.format(title=title))
                    
        if args.print_labels:
            output_file.write('# ------------------------------------------------------------------------------\n')
            wemdtools.bins.print_labels(region_set, output_file)
            output_file.write('# ------------------------------------------------------------------------------\n')

# Calculate averages and confidence intervals
sys.stdout.write('calculating averages and confidence intervals\n')
lifetimes = numpy.zeros((n_bins,6), numpy.float64)
eds = numpy.zeros((n_bins,n_bins,6), numpy.float64)
fluxes = numpy.zeros((n_bins,n_bins,6), numpy.float64)
rates  = numpy.zeros((n_bins,n_bins,6), numpy.float64)
transdat_ds = output_group['transitions']
transdat_iregion = transdat_ds['init_region'][:]
transdat_fregion = transdat_ds['final_region'][:]

def weighted_mean(tdata, field):
    norm = tdata['weight'].sum()
    wsum = (tdata[field] * tdata['weight']).sum()
    return wsum/norm

alpha = 1.0-args.confidence
n_sets = args.bssize or wemdtools.stats.mcbs.get_bssize(alpha)
syn_avg_eds = numpy.empty((n_sets,), numpy.float64)
syn_fluxes  = numpy.empty((n_sets,), numpy.float64)
syn_rates   = numpy.empty((n_sets,), numpy.float64)
lbi = int(math.floor(n_sets*alpha/2))
ubi = int(math.ceil(n_sets*(1-alpha/2)))
for iregion in xrange(0, n_bins):
    trans_iregion = transdat_ds[transdat_iregion == iregion]
    trans_has_lifetime = trans_iregion[trans_iregion['lifetime'] > 0]
    trans_has_duration = trans_iregion[trans_iregion['ed'] > 1]
    
    if not len(trans_iregion): continue
    
    # Get confidence limits on lifetime    
    lifetimes[iregion] = bootstrap_ci(weighted_mean, trans_has_lifetime[['lifetime', 'weight']], 
                                      alpha, n_sets, extended_output=True, args=('lifetime',))

    for fregion in xrange(0, n_bins):
        trans_ifregions = trans_has_duration[trans_has_duration['final_region'] == fregion]
        
        # Filter out data points where initial region probability is zero
        # This can happen, rarely, for fluctuations into and out of a target state 
        # (or it's an off-by-one error in the analysis somewhere...)
        trans_ifregions = trans_ifregions[trans_ifregions['iprob'] > 0]
        
        if not len(trans_ifregions): 
            continue
        
        this_total_weight         = trans_ifregions['weight'].sum()        
        fluxes[iregion,fregion,0] = this_total_weight
        rates[iregion,fregion,0]  = (trans_ifregions['weight'] / trans_ifregions['iprob']).sum()
        eds[iregion,fregion,0]    = (trans_ifregions['ed'] * trans_ifregions['weight'] ).sum() / this_total_weight
        
        # Use the same synthetic data sets for ed, flux, and rate CI calculations
        N = len(trans_ifregions)
        for iset in xrange(0, n_sets):
            indices = numpy.random.randint(N, size=(N,))
            syn_weights = trans_ifregions['weight'][indices]
            syn_iprob   = trans_ifregions['iprob'][indices]
            syn_eds     = trans_ifregions['ed'][indices]
            
            tot_weight = syn_weights.sum()
            
            syn_avg_eds[iset] = (syn_eds*syn_weights).sum() / tot_weight
            syn_fluxes[iset]  = tot_weight
            syn_rates[iset]   = (syn_weights / syn_iprob).sum()
            
        syn_avg_eds.sort()
        syn_fluxes.sort()
        syn_rates.sort()
        
        eds[iregion,fregion,1] = syn_avg_eds[lbi]; eds[iregion,fregion,2] = syn_avg_eds[ubi]
        fluxes[iregion,fregion,1] = syn_fluxes[lbi]; fluxes[iregion,fregion,2] = syn_fluxes[ubi]
        rates[iregion,fregion,1] = syn_rates[lbi]; rates[iregion,fregion,2] = syn_rates[ubi]
        
        for arr in fluxes, rates, eds:
            val, lb, ub = arr[iregion,fregion,0:3]
            arr[iregion,fregion,3] = ub - lb                            # CI width
            arr[iregion,fregion,4] = abs(arr[iregion,fregion,3] / val)  # relative width
            arr[iregion,fregion,5] = 2*max(ub-val,val-lb)               # symmetrized width

# Divide fluxes and rates by tau and number of tau considered, since we sum
# over start_iter to stop_iter; skip the relative width, with which dt divides out anyway
fluxes[...,[0,1,2,3,5]] /= args.dt * (stop_iter-start_iter+1) * (pcoord_len-1)
rates[...,[0,1,2,3,5]]  /= args.dt * (stop_iter-start_iter+1) * (pcoord_len-1)

# Multiply times else by dt
lifetimes[...,[0,1,2,3,5]] *= args.dt
eds[...,[0,1,2,3,5]] *= args.dt

# Write output files
max_ibin_width = len(str(n_bins-1))
# In order: initial bin, final bin, mean, lb, ub, width, rwidth, swidth
format_2d = '{ibin:{mw}d}    {fbin:{mw}d}    {0:20.15g}    {1:20.15g}    {2:20.15g}    {3:20.15g}    {4:20.15g}    {5:20.15g}\n'
format_1d = '{ibin:{mw}d}    {0:20.15g}    {1:20.15g}    {2:20.15g}    {3:20.15g}    {4:20.15g}    {5:20.15g}\n'
for iregion in xrange(0, n_bins):
    args.lifetime_stats.write(format_1d.format(*map(float,lifetimes[iregion]),ibin=iregion, mw=max_ibin_width))
    for fregion in xrange(0, n_bins):
        args.ed_stats.write(format_2d.format(*map(float,eds[iregion,fregion]),ibin=iregion,fbin=fregion, mw=max_ibin_width))
        args.flux_stats.write(format_2d.format(*map(float,fluxes[iregion,fregion]),ibin=iregion,fbin=fregion, mw=max_ibin_width))
        args.rate_stats.write(format_2d.format(*map(float,rates[iregion,fregion]),ibin=iregion,fbin=fregion, mw=max_ibin_width))

sys.stdout.write('\nFinal results:\n')
print('\nNumber of boundary crossings:',file=sys.stdout)
print_2d_int_array(transacc.n_crossings,file=sys.stdout)
print('Number of transitions:',file=sys.stdout)
print_2d_int_array(transacc.n_completions,file=sys.stdout)
