from __future__ import print_function, division
import os, sys, argparse, math
import numpy
import wemd
from wemd.util import extloader

import logging
log = logging.getLogger('w_ttimes')


def run_transanl_bf(args, sim_manager, region_set, pcoords):
    #from wemd.pcoords import RectilinearRegionSet
    #region_set = RectilinearRegionSet([[0, 0.4, 1.0, float('inf')]])
    region_set.clear()
    all_bins = region_set.get_all_bins()
    nbins = len(all_bins)
    all_indices = list(xrange(0,nbins))
    
    last_crossing   = numpy.zeros((nbins,nbins), numpy.uintp)
    last_entry      = numpy.zeros((nbins,), numpy.uintp)
    last_exit       = numpy.zeros((nbins,), numpy.uintp)
    last_completion = numpy.zeros((nbins,nbins), numpy.uintp)
    
    n_crossings     = numpy.zeros((nbins,nbins), numpy.uintp)
    n_completions   = numpy.zeros((nbins,nbins), numpy.uintp)
    
    dwell_weight_acc = numpy.zeros((nbins,), numpy.float64)
    dwell_acc       = numpy.zeros((nbins,), numpy.float64)
    dwell_sq_acc    = numpy.zeros((nbins,), numpy.float64)
    
    ed_weight_acc   = numpy.zeros((nbins,nbins), numpy.float64)
    ed_acc          = numpy.zeros((nbins,nbins), numpy.float64)
    ed_sq_acc       = numpy.zeros((nbins,nbins), numpy.float64)
    fpt_weight_acc  = numpy.zeros((nbins,nbins), numpy.float64)
    fpt_acc         = numpy.zeros((nbins,nbins), numpy.float64)
    fpt_sq_acc      = numpy.zeros((nbins,nbins), numpy.float64)

    
    tfile = args.trans_file
    edfile = args.ed_file 
    fptfile = args.fpt_file
    dwellfile = args.dwell_file
    
    if not args.suppress_headers:
        tfile.write('''\
# Crossings
# column 0: time point of end of transition
# column 1: initial->final bin indices
# column 2: weight at end of transition 
''')
        dwellfile.write('''\
# Dwell times (lifetimes)
# column 0: time point
# column 1: transition terminating dwell time
# column 2: bin index
# column 3: dwell time
# column 4: weight at end of transition
''')
        edfile.write('''\
# Event durations
# column 0: time point of end of transition
# column 1: initial->final bin indices
# column 2: event duration
# column 3: weight at end of transition
''')
        fptfile.write('''\
# First passage times
# column 0: time point of end of transition
# column 1: initial->final bin indices
# column 2: first passage time
# column 3: weight at end of transition
''')
        
    MIN=args.start
    MAX=args.stop or len(pcoords)
    CHUNK=args.chunksize
    output_file = args.output_file
    output_is_tty = output_file.isatty()
    
    last_index = region_set.map_to_indices([pcoords[MIN]])[0]
    for j in xrange(MIN+1, MAX, CHUNK):
        if output_is_tty:
            print('{:d}/{:d}'.format(j, MAX),file=output_file)
        pcoords_chunk = numpy.expand_dims(pcoords[j:j+CHUNK],1)
        indices_chunk = region_set.map_to_indices(pcoords_chunk)
        
        kmax = CHUNK if j+CHUNK <= MAX else MAX-j
        
        for k in xrange(0, kmax, 1):
            i = j+k
            weight = 1
            current_index = indices_chunk[k]
        
            if current_index != last_index:
                n_crossings[last_index,current_index] += 1
                
                clabel = '{:d}->{:d}'.format(long(last_index), long(current_index))
                
                tfile.write('{time:20d}    {label:<24s}    {weight:20.14e}\n'.format(time=i, label=clabel, weight=weight))
                
                if last_entry[last_index] > 0:
                    dwell = i - last_entry[last_index]
                    fdwell = float(dwell)
                    dwell_weight_acc[last_index] += weight
                    dwell_acc[last_index] += fdwell*weight
                    dwell_sq_acc[last_index] += fdwell*fdwell*weight 
                    dwellfile.write('{time:20d}    {label:<24s}    {dlabel:10s}    {dwell:20d}    {weight:20.14e}\n'
                                    .format(time=i, label=clabel, dlabel=str(long(last_index)), dwell=long(dwell), weight=weight))
                        
                # event duration X->Y: (time of this entry into Y) - (time of last exit from X)
                for istart in all_indices:
                    if istart == current_index: continue    # renewal
                    if istart == last_index: continue       # just a crossing
                    tlabel = '{:d}->{:d}'.format(long(istart), long(current_index))
                    
                    # initial region has been visited since last completion of initial->current;
                    # this indicates a transition from initial->current
                    if last_entry[istart] > last_completion[istart, current_index]:
                        # We can always compute event durations for nonadjacent states
                        ed = i - last_exit[istart]
                        fed = float(ed)
                        ed_weight_acc[istart,current_index] += weight
                        ed_acc[istart,current_index] += fed*weight
                        ed_sq_acc[istart,current_index] += fed*fed*weight
                        edfile.write('{time:20d}    {label:<24s}    {ed:20d}    {weight:20.14e}\n'
                                     .format(time=i, label=tlabel, ed=long(ed), weight=weight))
                        
                        # If we have seen a current->initial transition, then we can compute an FPT
                        # for initial->current
                        if last_completion[current_index, istart] > 0:
                            fpt = i - last_completion[current_index, istart]
                            ffpt = float(fpt)
                            fpt_weight_acc[istart,current_index] += weight
                            fpt_acc[istart,current_index] += ffpt*weight
                            fpt_sq_acc[istart,current_index] += ffpt*ffpt*weight
                            fptfile.write('{time:20d}    {label:20s}    {fpt:20d}    {weight:20.14e}\n'
                                         .format(time=i, label=tlabel, fpt=long(fpt), weight=weight))
                            
                        last_completion[istart,current_index] = i
                        n_completions[istart,current_index] += 1
                    
                last_exit[last_index] = i
                last_entry[current_index] = i
                last_crossing[last_index, current_index] = i
                    
            last_index = current_index
    
    
    avg_dwell = dwell_acc / dwell_weight_acc
    stdev_dwell = (dwell_sq_acc/dwell_weight_acc - avg_dwell*avg_dwell)**0.5
    avg_ed = ed_acc / ed_weight_acc
    stdev_ed = (ed_sq_acc/ed_weight_acc - avg_ed*avg_ed)**0.5
    avg_fpt = fpt_acc / fpt_weight_acc
    stdev_fpt = (fpt_sq_acc/fpt_weight_acc - avg_fpt*avg_fpt)**0.5
    
    print('Number of crossings:',file=output_file)
    print(n_crossings,file=output_file)
    print('Number of completed transitions:',file=output_file)
    print(n_completions,file=output_file)
    
    print('Average dwell times:',file=output_file)
    print(avg_dwell,file=output_file)
    print('Standard deviation of dwell times:',file=output_file)
    print(stdev_dwell,file=output_file)
    
    print('Average event durations:',file=output_file)
    print(avg_ed,file=output_file)
    print('Standard deviation of event durations:',file=output_file)
    print(stdev_ed,file=output_file)
    
    print('Average first passage times:',file=output_file)
    print(avg_fpt,file=output_file)
    print('Standard deviation of first passage times:',file=output_file)
    print(stdev_fpt,file=output_file)



parser = wemd.rc.common_arg_parser('w_ttimes', description='''
Perform lifetime and transition analysis on WEMD or brute force data.''')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
#parser.add_argument('-I', '--save-converted-input', )
parser.add_argument('-T', '--transitions', dest='trans_file', type=argparse.FileType('wt'), default='transitions.txt',
                    help='List transitions to TRANS_FILE (default: transitions.txt)')
parser.add_argument('-E', '--durations', dest='ed_file', type=argparse.FileType('wt'), default='durations.txt',
                    help='List event durations to ED_FILE (default: durations.txt)')
parser.add_argument('-F', '--fpts', dest='fpt_file', type=argparse.FileType('wt'), default='fpts.txt',
                    help='List first passage times to FPT_FILE (default: fpts.txt)')
parser.add_argument('-D', '--dwells', dest='dwell_file', type=argparse.FileType('wt'), default='dwells.txt',
                    help='List dwell times to DWELL_FILE (default: dwells.txt)')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
parser.add_argument('-b', '--bins', dest='binexpr',
                    help='''Use BINEXPR to construct bins for analysis; must contain an assignment to 
                     a variable called "region_set" (Default: load from system file.)''')
parser.add_argument('--start', dest='start', type=int, default=0,
                    help='Start at entry START (zero-based, default: 0)')
parser.add_argument('--stop', dest='stop', type=int,
                    help='Stop before entry STOP (zero-based, default: process all data)')
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
    runtime_config = wemd.rc.read_config(args.run_config_file)
    runtime_config.update_from_object(args)
    sim_manager = wemd.rc.load_sim_manager(runtime_config)
    sim_manager.load_data_manager()
    sim_manager.data_manager.open_backing()
    sim_manager.load_system_driver()
    region_set = sim_manager.system.region_set
    
if args.datafile:
    print('reading data from', args.datafile, file=args.output_file)
    pcoords = numpy.load(args.datafile, 'r')
    wemd.rc.default_cmdline_dispatch(run_transanl_bf, 
                                     args=[args, sim_manager, region_set, pcoords], kwargs={}, cmdline_args=args, log=log)
else:
    print('reading data from WEMD simulation as specified in', args.run_config_file, file=args.output_file)
    