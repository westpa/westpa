from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse
import numpy
import wemd, wemdtools
from wemdtools.transitions.transacc import TransitionEventAccumulator
from wemdtools.trajectories.trajtree import TrajTree 
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
            args.output_file.write('{:>{nsegwidth}d} of a maximum of {:<{nsegwidth}d} segments analyzed\n'
                                   .format(n_visited, n_to_visit, nsegwidth=nsegwidth))
        else:
            args.output_file.write('{:>{nsegwidth}d} of {:<{nsegwidth}d} segments ({:5.1f}%) analyzed\n'
                                   .format(n_visited, n_to_visit, pct_visited, nsegwidth=nsegwidth))
        args.output_file.flush()

    if len(history) == 0:
        transacc.accumulate_transitions(segment.pcoord, weight=segment.weight, 
                                        region_weights=binprobs[segment.n_iter-start_iter],
                                        continuation = False, label=str(segment.n_iter))
    else:
        transacc.accumulate_transitions(segment.pcoord[1:], weight=segment.weight, 
                                        region_weights=binprobs[segment.n_iter-start_iter], 
                                        continuation = True, label=str(segment.n_iter))
        
    # Delete the progress coordinate from processed segments so that memory use doesn't explode
    del segment.pcoord
    segment.pcoord = None


parser = wemd.rc.common_arg_parser('w_ttimes', description='''
Perform lifetime, transition, and kinetic analysis on WEMD data.''')
parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                    help='Do not emit status messages (default: emit status every 1000 WE segments)')

parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)

parser.add_argument('-T', '--transitions', dest='trans_file', type=argparse.FileType('wt'),
                    help='List transitions to TRANS_FILE (default: do not list)')
parser.add_argument('-E', '--eds', dest='ed_file', type=argparse.FileType('wt'), default='eds.txt',
                    help='List event durations to ED_FILE (default: durations.txt)')
parser.add_argument('-R', '--rates', dest='rate_file', type=argparse.FileType('wt'), default='rates.txt',
                    help='List bin-to-bin rates to RATE_FILE (default: rates.txt)')
parser.add_argument('-F', '--fpts', dest='fpt_file', type=argparse.FileType('wt'),
                    help='List first passage times to FPT_FILE (default: not analyzed)')
parser.add_argument('-L', '--lifetimes', dest='lifetime_file', type=argparse.FileType('wt'),
                    help='List times spent in each region to LIFETIME_FILE (default: do not list)')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')

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

parser.add_argument('datafile', nargs='*',
                    help='Read progress coordinate from DATAFILE(s) (default: load WEMD HDF5 file specified in wemd.cfg).')
args = parser.parse_args()

wemd.rc.config_logging(args, 'w_ttimes')

runtime_config = None
sim_manager = None
region_set = wemdtools.bins.get_region_set_from_args(args, status_stream=sys.stderr)
n_bins = len(region_set.get_all_bins())

if not args.suppress_headers:
    if args.trans_file is not None:
        args.trans_file.write('''\
# Transitions
# column 0: label (trajectory for brute force data, iteration number for WE data)
# column 1: time point of end of transition
# column 2: initial->final bin indices
# column 3: weight at end of transition
# (column widths:  [24, 20, 24, 20])
# (column offsets assuming no field overflows: [0, 28, 52, 80]) 
''')
    if args.lifetime_file is not None:
        args.lifetime_file.write('''\
# Lifetimes
# column 0: label (trajectory for brute force data, iteration number for WE data)
# column 1: time point
# column 2: bin index
# column 3: lifetime
# column 4: weight at end of transition
# (column widths:  [24, 20, 10, 20, 20])
# (column offsets assuming no field overflows: [0, 28, 52, 66, 90] ) 
''')
    if args.ed_file is not None:
        args.ed_file.write('''\
# Event durations
# column 0: (trajectory for brute force data, iteration number for WE data)
# column 1: time point of end of transition
# column 2: initial->final bin indices
# column 3: event duration
# column 4: weight at end of transition
# (column widths:  [24, 20, 24, 20, 20])
# (column offsets assuming no field overflows: [0, 28, 52, 80, 104]) 
''')
    if args.fpt_file is not None:
        args.fpt_file.write('''\
# First passage times
# column 0: (trajectory for brute force data, iteration number for WE data)
# column 1: time point of end of transition
# column 2: initial->final bin indices
# column 3: first passage time
# column 4: weight at end of transition
# (column widths:  [24, 20, 24, 20, 20])
# (column offsets assuming no field overflows: [0, 28, 52, 80, 104]) 
''')
    if args.rate_file is not None:
        args.rate_file.write('''\
# Kinetic rates
# column 0: (trajectory for brute force data, iteration number for WE data)
# column 1: time point of end of transition
# column 2: initial->final bin indices
# column 3: observed rate
''')
    

if args.datafile:
    from wemd.util.config_dict import ConfigDict
    if runtime_config is None: runtime_config = ConfigDict()
    runtime_config['data.h5file'] = args.datafile[0]
    print('reading WEMD data from', args.datafile[0], file=args.output_file)
else:
    print('reading data from WEMD simulation', file=args.output_file)
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
    args.output_file.write('start and stop iterations are the same -- no data to analyze')
    sys.exit(0)

if args.whole_only:
    args.output_file.write('considering whole trajectories only\n')

if args.cache_pcoords:
    args.output_file.write('will cache pcoord data in memory\n')
    

tree = TrajTree(sim_manager.data_manager, cache_pcoords = args.cache_pcoords)
transacc = TransitionEventAccumulator(region_set, args.trans_file, args.lifetime_file, args.ed_file, args.fpt_file, args.rate_file)

# Use the caching wrapper, which is conveniently created around the true data manager when instantiating TrajTree
data_manager = tree.data_manager

# Determine the shape of pcoords - (n_segs, pcoord_len, pcoord_ndim)
pcoord_ds = data_manager.get_pcoord_dataset(1)
pcoord_len = pcoord_ds.shape[1]
pcoord_ndim = pcoord_ds.shape[2]

if not args.quiet_mode:
    args.output_file.write('determining bin probabilities\n')
binprobs = numpy.empty((stop_iter-start_iter+1, pcoord_len, n_bins),dtype=numpy.float64)
for n_iter in xrange(start_iter, stop_iter+1):
    if not args.quiet_mode and (n_iter - start_iter) % 10 == 0:
        args.output_file.write('iteration {}\n'.format(n_iter))
    i_iter = n_iter - start_iter
    pcoords = data_manager.get_pcoord_array(n_iter)
    weights = data_manager.get_seg_index(n_iter)['weight']
    for ti in xrange(0,pcoord_len):
        index_map = region_set.map_to_all_indices(pcoords[:,ti,:])
        for si in xrange(0, len(index_map)):
            binprobs[i_iter,ti,index_map[si]] += weights[si]
        assert abs(binprobs[i_iter,ti,:].sum() - 1.0) < 1.0e-8    

args.output_file.write('tracing trajectories\n')
tree.trace_trajectories(start_iter, stop_iter, 
                        callable=accumulate_transitions_segment,
                        get_state=transacc.get_state,
                        set_state=transacc.set_state,
                        args=(transacc,data_manager,tree,binprobs,start_iter,args),
                        whole_only = bool(args.whole_only))
args.output_file.write('{} segments analyzed\n'.format(tree.segments_visited))

args.output_file.write('Final results:\n')
print('\nNumber of crossings:',file=args.output_file)
print_2d_int_array(transacc.n_crossings,file=args.output_file)
print('Number of completed long-time/non-adjacent transitions:',file=args.output_file)
print_2d_int_array(transacc.n_completions,file=args.output_file)

if transacc.track_lifetimes:
    print('\nAverage lifetimes:',file=args.output_file)
    print_2d_float_array(numpy.expand_dims(transacc.lt_acc.average(),1),file=args.output_file)
    print('Standard deviation of lifetimes:',file=args.output_file)
    print_2d_float_array(numpy.expand_dims(transacc.lt_acc.std(),1),file=args.output_file)

if transacc.track_eds:
    print('\nAverage event durations:',file=args.output_file)
    print_2d_float_array(transacc.ed_acc.average(),file=args.output_file)
    print('Standard deviation of event durations:',file=args.output_file)
    print_2d_float_array(transacc.ed_acc.std(),file=args.output_file)
    
if transacc.track_rates:
    print('\nAverage fluxes:', file=args.output_file)
    print_2d_float_array(transacc.flux_acc.average(), file=args.output_file)
    print('Standard deviation of fluxes:', file=args.output_file)
    print_2d_float_array(transacc.flux_acc.std(), file=args.output_file)
    
    print('\nAverage rates:', file=args.output_file)
    print_2d_float_array(transacc.rate_acc.average(), file=args.output_file)
    print('Standard deviation of rates:',file=args.output_file)
    print_2d_float_array(transacc.rate_acc.std(),file=args.output_file)

if transacc.track_fpts:    
    print('\nAverage first passage times:',file=args.output_file)
    print_2d_float_array(transacc.fpt_acc.average(),file=args.output_file)
    print('Standard deviation of first passage times:',file=args.output_file)
    print_2d_float_array(transacc.fpt_acc.std(),file=args.output_file)
