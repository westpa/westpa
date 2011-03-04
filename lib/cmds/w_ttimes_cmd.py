from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math
from collections import namedtuple
from itertools import count, izip, islice
import numpy
import wemd
from wemd.util import extloader
from wemdtools.transitions.transacc import TransitionEventAccumulator 
from math import ceil, log10

import logging
log = logging.getLogger('w_ttimes')

def print_2d_int_array(array, format='{:6d}', sep='  ', file=sys.stdout):
    for row in array:
        file.write('{}\n'.format(sep.join(format.format(long(field)) for field in row)))
    try:
        file.flush()
    except AttributeError:
        pass

def print_2d_float_array(array, format='{:12.6e}', sep='  ', file=sys.stdout):
    for row in array:
        file.write('{}\n'.format(sep.join(format.format(float(field)) for field in row)))
    try:
        file.flush()
    except AttributeError:
        pass
        
            
def report_counts_averages(args, transacc, output_file):
    print('Number of crossings:',file=output_file)
    print_2d_int_array(transacc.n_crossings,file=output_file)
    print('Number of completed transitions:',file=output_file)
    print_2d_int_array(transacc.n_completions,file=output_file)
    
    if transacc.track_lifetimes:
        print('Average lifetimes:',file=output_file)
        print_2d_float_array(numpy.expand_dims(transacc.lt_acc.average(),1),file=output_file)
        print('Standard deviation of dwell times:',file=output_file)
        print_2d_float_array(numpy.expand_dims(transacc.lt_acc.std(),1),file=output_file)
    
    if transacc.track_eds:
        print('Average event durations:',file=output_file)
        print_2d_float_array(transacc.ed_acc.average(),file=output_file)
        print('Standard deviation of event durations:',file=output_file)
        print_2d_float_array(transacc.ed_acc.std(),file=output_file)

    if transacc.track_fpts:    
        print('Average first passage times:',file=output_file)
        print_2d_float_array(transacc.fpt_acc.average(),file=output_file)
        print('Standard deviation of first passage times:',file=output_file)
        print_2d_float_array(transacc.fpt_acc.std(),file=output_file)
    
    
def run_transanl_bf(args, transacc, pcoords):
    maxidx = args.stop or len(pcoords)
    n_chunks = int(math.ceil(maxidx/args.chunksize))
    for ichunk in xrange(0, n_chunks):
        lbi = ichunk*args.chunksize
        ubi = min((ichunk+1)*args.chunksize, maxidx)
        pc_chunk = pcoords[lbi:ubi]
        if pc_chunk.ndim == 1:
            pc_chunk = numpy.expand_dims(pc_chunk,-1)

        if not args.quiet_mode and ichunk > 0:    
            args.output_file.write('{current_time:d}/{total_time:d}\n'.format(current_time=lbi,
                                                                              total_time=maxidx))
            report_counts_averages(args, transacc, args.output_file)
            args.output_file.write('\n')
            args.output_file.flush()

        if ichunk == 0:
            transacc.accumulate_transitions(pc_chunk, continuation=False)
        else:
            transacc.accumulate_transitions(pc_chunk, continuation=True)

def accumulate_transitions_segment(segment, children, history, transacc, data_manager, trajtree, args):
    n_visited = trajtree.segments_visited
    n_to_visit = trajtree.segments_to_visit
    if not args.quiet_mode and n_visited % 1000 == 0:
        nsegwidth = int(ceil(log10(n_to_visit)))
        pct_visited = n_visited / n_to_visit * 100
        args.output_file.write('{:>{nsegwidth}d} of {:<{nsegwidth}d} segments ({:5.1f}%) analyzed\n'
                               .format(n_visited, n_to_visit, pct_visited, nsegwidth=nsegwidth))
        args.output_file.flush()
    if segment.p_parent_id < 0:
        transacc.accumulate_transitions(segment.pcoord, weight=segment.weight, continuation = False)
    else:
        transacc.accumulate_transitions(segment.pcoord[1:], weight=segment.weight, continuation = True)
    
def run_transanl_wemd(args, transacc, sim_manager):
    start_iter = args.start or 1
    if not args.stop: args.stop = sim_manager.data_manager.current_iteration - 1
    stop_iter  = min(args.stop, sim_manager.data_manager.current_iteration-1)
    if start_iter == stop_iter: return
    
    from wemdtools.trajectories.trajtree import TrajTree
    tree = TrajTree(sim_manager.data_manager)
    tree.trace_trajectories(start_iter, stop_iter, 
                            callable=accumulate_transitions_segment,
                            get_state=transacc.get_state,
                            set_state=transacc.set_state,
                            args=(transacc,sim_manager.data_manager,tree,args))
    args.output_file.write('{} segments analyzed\n'.format(tree.segments_visited))
    

parser = wemd.rc.common_arg_parser('w_ttimes', description='''
Perform lifetime and transition analysis on WEMD or brute force data.''')
parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                    help='Do not emit intermediate results (default: emit intermediate results '
                        +'every CHUNKSIZE time points)')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
parser.add_argument('-T', '--transitions', dest='trans_file', type=argparse.FileType('wt'),
                    help='List transitions to TRANS_FILE (default: do not list)')
parser.add_argument('-E', '--eds', dest='ed_file', type=argparse.FileType('wt'), default='eds.txt',
                    help='List event durations to ED_FILE (default: durations.txt)')
parser.add_argument('-F', '--fpts', dest='fpt_file', type=argparse.FileType('wt'),
                    help='List first passage times to FPT_FILE (default: not analyzed)')
parser.add_argument('-L', '--lifetimes', dest='lifetime_file', type=argparse.FileType('wt'),
                    help='List times spent in each region to LIFETIME_FILE (default: do not list)')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
parser.add_argument('-b', '--bins', dest='binexpr',
                    help='''Use BINEXPR to construct bins for analysis; must contain an assignment to 
                     a variable called "region_set" (Default: load from system file.)''')
parser.add_argument('--start', dest='start', type=int, default=0,
                    help='Start at brute force entry or WEMD iteration START; ' 
                        +'zero-based for brute force, one-based for WEMD (default: the beginning)')
parser.add_argument('--stop', dest='stop', type=int,
                    help='Stop at brute force entry or WEMD iteration STOP; '
                        +'zero-based for brute force, one-based for WEMD (default: process all data)')
parser.add_argument('-C', '--chunksize', dest='chunksize', type=int, default=100000,
                    help='Retrieve and analyze data in chunks of size CHUNKSIZE (default: 100000)')
parser.add_argument('datafile', nargs='?',
                    help='Read progress coordinate from DATAFILE (default: load from WEMD HDF5 file).')
args = parser.parse_args()

wemd.rc.config_logging(args, 'w_ttimes')

runtime_config = None
sim_manager = None
if args.binexpr:
    print('using bins specified on the command line',file=args.output_file)
    namespace = globals()
    namespace.update(locals())
    try:
        ccode = compile(args.binexpr, '<string>', mode='exec')
        exec ccode in namespace
        region_set = namespace['region_set']
    except SyntaxError as e:
        sys.stderr.write('{!s}\n'.format(e))
        sys.exit(1)
    except KeyError:
        sys.stderr.write('inline region set code MUST assign to a variable named "region_set"\n')
        sys.exit(1)
else:
    print('loading bin data from WEMD system', file=args.output_file)
    runtime_config = wemd.rc.read_config(args.run_config_file)
    runtime_config.update_from_object(args)
    sim_manager = wemd.rc.load_sim_manager(runtime_config)
    sim_manager.load_system_driver()
    region_set = sim_manager.system.region_set

if not args.suppress_headers:
    if args.trans_file is not None:
        args.trans_file.write('''\
# Crossings
# column 0: time point of end of transition
# column 1: initial->final bin indices
# column 2: weight at end of transition 
''')
    if args.lifetime_file is not None:
        args.lifetime_file.write('''\
# Lifetimes
# column 0: time point
# column 1: bin index
# column 2: dwell time
# column 3: weight at end of transition
''')
    if args.ed_file is not None:
        args.ed_file.write('''\
# Event durations
# column 0: time point of end of transition
# column 1: initial->final bin indices
# column 2: event duration
# column 3: weight at end of transition
''')
    if args.fpt_file is not None:
        args.fpt_file.write('''\
# First passage times
# column 0: time point of end of transition
# column 1: initial->final bin indices
# column 2: first passage time
# column 3: weight at end of transition
''')
    
transacc = TransitionEventAccumulator(region_set, args.trans_file, args.lifetime_file, args.ed_file, args.fpt_file)

if args.datafile:
    print('reading data from', args.datafile, file=args.output_file)
    pcoords = numpy.load(args.datafile, 'r')
    wemd.rc.default_cmdline_dispatch(run_transanl_bf, 
                                     args=[args, transacc, pcoords], kwargs={}, cmdline_args=args, log=log)
else:
    print('reading data from WEMD simulation', file=args.output_file)
    if runtime_config is None:
        runtime_config = wemd.rc.read_config(args.run_config_file)
        runtime_config.update_from_object(args)        
    if sim_manager is None:
        sim_manager = wemd.rc.load_sim_manager(runtime_config)
    sim_manager.load_data_manager()
    sim_manager.data_manager.open_backing(mode='r')
    wemd.rc.default_cmdline_dispatch(run_transanl_wemd, args=[args, transacc, sim_manager], kwargs={}, cmdline_args=args, log=log)
    
args.output_file.write('Final results:\n')
report_counts_averages(args, transacc, args.output_file)
    