from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math
from collections import namedtuple
from itertools import count, izip, islice
import numpy
import wemd
from wemd.pcoords import PiecewiseRegionSet, RectilinearRegionSet
from wemd.util import extloader
from wemdtools.transitions.transacc import TransitionEventAccumulator
from wemdtools.files import load_npy_or_text 
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
        
            
def report_counts_averages(args, transacc, output_file):
    print('Number of crossings:',file=output_file)
    print_2d_int_array(transacc.n_crossings,file=output_file)
    print('Number of completed transitions:',file=output_file)
    print_2d_int_array(transacc.n_completions,file=output_file)
    
    if transacc.track_lifetimes:
        print('Average lifetimes:',file=output_file)
        print_2d_float_array(numpy.expand_dims(transacc.lt_acc.average(),1),file=output_file)
        print('Standard deviation of lifetimes:',file=output_file)
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
    
    
def run_transanl_bf(args, transacc):
    for (ifile, datafile) in enumerate(args.datafile):
        print('reading data from', datafile, file=args.output_file)
        pcoords = load_npy_or_text(datafile)
        if args.ignore_1st_col:
            pcoords = pcoords[:,1:]
    
        minidx = args.start or 0
        maxidx = args.stop or len(pcoords)
        n_chunks = int(math.ceil((maxidx-minidx)/args.chunksize))
        for ichunk in xrange(0, n_chunks):
            lbi = minidx+ichunk*args.chunksize
            ubi = min(minidx+(ichunk+1)*args.chunksize, maxidx)
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
                transacc.accumulate_transitions(pc_chunk, continuation=False, label=str(ifile))
            else:
                transacc.accumulate_transitions(pc_chunk, continuation=True, label=str(ifile))
        del pcoords

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
        transacc.accumulate_transitions(segment.pcoord, weight=segment.weight, continuation = False, label=str(segment.n_iter))
    else:
        transacc.accumulate_transitions(segment.pcoord[1:], weight=segment.weight, continuation = True, label=str(segment.n_iter))
        
    # Delete the progress coordinate from processed segments so that memory use doesn't explode
    del segment.pcoord
    segment.pcoord = None
    
def run_transanl_wemd(args, transacc, sim_manager):
    start_iter = args.start or 1
    if not args.stop: args.stop = sim_manager.data_manager.current_iteration - 1
    stop_iter  = min(args.stop, sim_manager.data_manager.current_iteration-1)
    if start_iter == stop_iter: return
    
    from wemdtools.trajectories.trajtree import TrajTree
    tree = TrajTree(sim_manager.data_manager, cache_pcoords = args.cache_pcoords)
    
    if args.cache_pcoords:
        args.output_file.write('will cache pcoord data in memory\n')
    
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
                        +'every CHUNKSIZE brute force time points or every 1000 WE segments)')

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

parser.add_argument('--binexpr', dest='binexpr',
                    help='''Use BINEXPR to construct bins for analysis; must contain an assignment to 
                     a variable called "region_set" (Default: load from system file.)''')
parser.add_argument('--binbounds', dest='binbounds',
                    help='''Construct rectilinear bins from BINBOUNDS, which will be parsed as a list of lists
                    of bin boundaries.  The text 'inf' will be replaced by infinity.''')

parser.add_argument('--start', dest='start', type=int, default=0,
                    help='Start at brute force entry or WEMD iteration START; ' 
                        +'zero-based for brute force, one-based for WEMD (default: the beginning)')
parser.add_argument('--stop', dest='stop', type=int,
                    help='Stop at brute force entry or WEMD iteration STOP; '
                        +'zero-based for brute force, one-based for WEMD (default: process all data)')

parser.add_argument('-C', '--chunksize', dest='chunksize', type=int, default=100000,
                    help='Retrieve and analyze brute force data in chunks of size CHUNKSIZE (default: 100000)')
parser.add_argument('--ignore-first-column', dest='ignore_1st_col', action='store_true',
                    help='Ignore the first column in data files (useful for ignoring a time column). '
                        +'(Default: do not ignore.)')

parser.add_argument('--wemd', dest='force_wemd', action='store_true',
                    help = 'force interpretation of the given DATAFILE as a WEMD HDF5 file '
                         + '(so one doesn\'t have to create a wemd.cfg just to run this program)')
parser.add_argument('--cache-pcoords', dest='cache_pcoords', action='store_true',
                    help = '(WE only) cache progress coordinate data in memory for improved analysis speed; '
                         + 'generally requires as much free RAM as the size of the HDF5 file (default: do not cache')

parser.add_argument('datafile', nargs='*',
                    help='Read progress coordinate from DATAFILE(s) (default: load from WEMD HDF5 file). '
                        +'Multiple brute-force trajectories are supported. WEMD HDF5 files are supported '
                        +'by naming them AND also specifying the --wemd switch.')
args = parser.parse_args()

wemd.rc.config_logging(args, 'w_ttimes')

runtime_config = None
sim_manager = None
if args.binexpr:
    print('using region set expression on the command line',file=args.output_file)
    print('region set expression:',file=args.output_file)
    print(args.binexpr,file=args.output_file)
    print('----', file=args.output_file)
    namespace = globals()
    namespace.update(locals())
    try:
        ccode = compile(args.binexpr, '<BINEXPR>', mode='exec')
        exec ccode in namespace
        region_set = namespace['region_set']
    except SyntaxError as e:
        sys.stderr.write('{!s}\n'.format(e))
        sys.exit(1)
    except KeyError:
        sys.stderr.write('inline region set code MUST assign to a variable named "region_set"\n')
        sys.exit(1)
elif args.binbounds:
    import re
    re_transform_inf = re.compile(r'inf', re.IGNORECASE)
    
    transformed_expr = '{}\n'.format(re_transform_inf.sub("float('inf')", args.binbounds))
    binbounds = eval(transformed_expr)
    
    print('using rectiliear bin boundaries specified on the command line')
    print('bin expression:',file=args.output_file)
    print(binbounds,file=args.output_file)
    print('----', file=args.output_file)
    region_set = RectilinearRegionSet(binbounds) 
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
    
transacc = TransitionEventAccumulator(region_set, args.trans_file, args.lifetime_file, args.ed_file, args.fpt_file)

if args.datafile and not args.force_wemd:
    wemd.rc.default_cmdline_dispatch(run_transanl_bf, 
                                     args=[args, transacc], kwargs={}, cmdline_args=args, log=log)
else:
    if args.force_wemd:
        from wemd.util.config_dict import ConfigDict
        if runtime_config is None: runtime_config = ConfigDict()
        runtime_config['data.h5file'] = args.datafile
        print('reading WEMD data from', args.datafile, file=args.output_file)
    else:
        print('reading data from WEMD simulation', file=args.output_file)
        if runtime_config is None:
            runtime_config = wemd.rc.read_config(args.run_config_file)
                    
    if sim_manager is None:
        sim_manager = wemd.rc.load_sim_manager(runtime_config)
    sim_manager.load_data_manager()
    sim_manager.data_manager.open_backing(mode='r')
    wemd.rc.default_cmdline_dispatch(run_transanl_wemd, args=[args, transacc, sim_manager], kwargs={}, cmdline_args=args, log=log)
    
args.output_file.write('Final results:\n')
report_counts_averages(args, transacc, args.output_file)
    