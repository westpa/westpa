from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math
from collections import namedtuple
from itertools import count, izip, islice
import numpy
import wemd
from wemd.util import extloader

import logging
log = logging.getLogger('w_ttimes')

NAN=float('nan')

class RunningStatsAccumulator:
    def __init__(self, shape, dtype=numpy.float64, count_dtype=numpy.uint, weight_dtype=numpy.float64):
        self.sum = numpy.zeros(shape, dtype)
        self.sqsum = numpy.zeros(shape, dtype)
        self.weight = numpy.zeros(shape, weight_dtype)
        self.count  = numpy.zeros(shape, count_dtype)
        
    def incorporate(self, index, value, weight):
        self.count[index] += 1
        self.weight[index] += weight
        self.sum[index] += value
        self.sqsum[index] += value*value
        
    def average(self):
        valid = (self.weight > 0.0)
        avg = numpy.zeros_like(self.sum)
        avg[~valid] = NAN
        avg[valid] = self.sum[valid] / self.weight[valid]
        return avg
    mean = average
    
    def std(self):
        valid = (self.weight > 0.0)
        avg = self.average()
        sqavg = numpy.zeros_like(self.sqsum)
        sqavg[~valid] = NAN
        sqavg[valid]  = self.sqsum[valid] / self.weight[valid]
        return (sqavg - avg*avg)**0.5

class TransitionEventAccumulator:
    index_dtype  = numpy.uintp
    count_dtype  = numpy.uint64
    time_dtype   = numpy.float64
    weight_dtype = numpy.float64
    
    def __init__(self, region_set, tfile = None, ltfile = None, edfile = None, fptfile = None):
        self.region_set = region_set
        self.tfile = tfile
        self.edfile = edfile
        self.fptfile = fptfile
        self.ltfile = ltfile
        
        self.track_transitions = (tfile is not None)
        self.track_lifetimes = (ltfile is not None)
        self.track_fpts = (fptfile is not None)
        self.track_eds = (edfile is not None)
        
        self.n_bins = len(self.region_set.get_all_bins())
        self.clear()
        
    def clear(self):
        self._clear_arrays()
        
        self.n_crossings     = numpy.zeros((self.n_bins,self.n_bins), self.count_dtype)
        self.n_completions   = numpy.zeros((self.n_bins,self.n_bins), self.count_dtype)
        
        self.lt_acc  = RunningStatsAccumulator(shape=(self.n_bins,), dtype=self.time_dtype, count_dtype=self.count_dtype,
                                               weight_dtype = self.weight_dtype)
        self.ed_acc  = RunningStatsAccumulator(shape=(self.n_bins,self.n_bins), dtype=self.time_dtype, count_dtype=self.count_dtype,
                                               weight_dtype = self.weight_dtype)
        self.fpt_acc = RunningStatsAccumulator(shape=(self.n_bins,self.n_bins), dtype=self.time_dtype, count_dtype=self.count_dtype,
                                               weight_dtype = self.weight_dtype)
        
        self.time_index = None
        self.last_region = None        
        
    def _clear_arrays(self):
        self.last_crossing   = numpy.zeros((self.n_bins,self.n_bins), self.index_dtype)
        self.last_entry      = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_exit       = numpy.zeros((self.n_bins,), self.index_dtype)
        self.last_completion = numpy.zeros((self.n_bins,self.n_bins), self.index_dtype)
    
    def record_transition_entry(self, time_index, initial_region, final_region, weight):
        if self.tfile:
            label = '{:d}->{:d}'.format(int(initial_region),int(final_region))
            self.tfile.write('{time_index:20d}    {label:<24s}    {weight:20.14e}\n'
                             .format(time_index=long(time_index), label=label, weight=float(weight)))
    
    def record_lifetime_entry(self, time_index, final_region, lifetime, weight):
        #self.lt_acc.incorporate(final_region, lifetime, weight)
            #self.count[index] += 1
            #self.weight[index] += weight
            #self.sum[index] += value
            #self.sqsum[index] += value*value
        lt_acc = self.lt_acc
        lt_acc.count[final_region] += 1
        lt_acc.weight[final_region] += weight
        lt_acc.sum[final_region] += lifetime
        lt_acc.sqsum[final_region] += lifetime*lifetime
        
        if self.ltfile:
            self.ltfile.write('{time_index:20d}    {label:10d}    {lifetime:20d}    {weight:20.14e}\n'
                              .format(time_index=long(time_index), label=long(final_region), lifetime=long(lifetime), 
                                      weight=float(weight)))
    
    def record_ed_entry(self, time_index, initial_region, final_region, ed, weight):
        #self.ed_acc.incorporate((initial_region,final_region), ed, weight)
        ed_acc = self.ed_acc
        ed_acc.count[initial_region,final_region] += 1
        ed_acc.weight[initial_region,final_region] += weight
        ed_acc.sum[initial_region,final_region] += ed
        ed_acc.sqsum[initial_region,final_region] += ed*ed
        if self.edfile:
            label = '{:d}->{:d}'.format(int(initial_region),int(final_region))
            self.edfile.write('{time_index:20d}    {label:<24s}    {ed:20d}    {weight:20.14e}\n'
                              .format(time_index=long(time_index), label=label, ed=long(ed), weight=float(weight)))
    
    def record_fpt_entry(self, time_index, initial_region, final_region, fpt, weight):
        #self.fpt_acc.incorporate((initial_region,final_region), fpt, weight)
        fpt_acc = self.fpt_acc
        fpt_acc.count[initial_region,final_region] += 1
        fpt_acc.weight[initial_region,final_region] += weight
        fpt_acc.sum[initial_region,final_region] += fpt
        fpt_acc.sqsum[initial_region,final_region] += fpt*fpt
        if self.fptfile:
            label = '{:d}->{:d}'.format(int(initial_region),int(final_region))
            self.fptfile.write('{time_index:20d}    {label:<24s}    {fpt:20d}    {weight:20.14e}\n'
                              .format(time_index=long(time_index), label=label, fpt=long(fpt), weight=float(weight)))
    
    def accumulate_transitions(self, pcoords, weights = None, time_index=None, continuation = False):
        """Assign the given progress coordinates to regions, and then determine transitions among regions.
        If continuation is True, ignore the first point (raises an error if accumulate_transitions() has not
        been called at least once already), otherwise use the first point as the initial point."""
        
        if weights is not None:
            if len(weights) != len(pcoords):
                raise ValueError('weights array shape [{!r}] is not compatible with pcoords shape [{!r}]'
                                 .format(weights.shape, pcoords.shape))
        else:
            # weights is None
            weights = numpy.ones((len(pcoords),))
            
        
        indices = self.region_set.map_to_indices(pcoords)
        
        # a table of time points, region at that time point, and (yes, believe it), the region at the prior time point
        # in that order
        TIMEPOINT = 0; CURRENT_REGION = 1; LAST_REGION = 2
        
        if continuation:
            tracking_table = numpy.empty((len(pcoords), 3), self.index_dtype)
            tracking_table[0,LAST_REGION] = self.last_region
            tracking_table[:,CURRENT_REGION] = indices
            if len(pcoords) > 1:
                tracking_table[1:,LAST_REGION] = indices[:-1]
            tracking_table[:,TIMEPOINT] = numpy.arange(self.time_index+1, self.time_index+len(pcoords)+1, dtype=self.index_dtype)
        else:
            # not a continuation
            self.clear()
            self.time_index = time_index = time_index or 0
            tracking_table = numpy.empty((len(pcoords)-1,3), self.index_dtype)
            tracking_table[:,LAST_REGION] = indices[:-1]
            tracking_table[:,CURRENT_REGION] = indices[1:]
            tracking_table[:,TIMEPOINT] = numpy.arange(self.time_index+1, self.time_index+len(pcoords), dtype=self.index_dtype)
        
        observed_at = tracking_table[:,CURRENT_REGION] != tracking_table[:,LAST_REGION]
        trans_observed = tracking_table[observed_at, :]
        weights_observed = weights[observed_at]

        # pull all references to self out of the inner loops for a 25% reduction in CPU time!
        last_crossing = self.last_crossing
        last_entry = self.last_entry
        last_exit = self.last_exit
        last_completion = self.last_completion
        n_crossings = self.n_crossings
        n_completions = self.n_completions
        record_transition_entry = self.record_transition_entry
        record_lifetime_entry = self.record_lifetime_entry
        record_ed_entry = self.record_ed_entry
        record_fpt_entry = self.record_fpt_entry
        n_bins = self.n_bins
        index_dtype = self.index_dtype
        track_transitions = self.track_transitions
        track_lifetimes = self.track_lifetimes
        track_eds = self.track_eds
        track_fpts = self.track_fpts

        irdisc = numpy.zeros((n_bins,), numpy.bool_)        
        init_regions = numpy.arange(0,n_bins,dtype=index_dtype)
        
        for ((time_index,current_region,last_region), weight) in izip(trans_observed, weights_observed):
            n_crossings[last_region,current_region] += 1
            if track_transitions:
                record_transition_entry(time_index,last_region,current_region,weight)
            
            if track_lifetimes and last_entry[last_region] > 0:
                lifetime = time_index - last_entry[last_region]
                record_lifetime_entry(time_index, last_region, lifetime, weight)
                
            irdisc[:] = True
            irdisc[last_region] = False
            irdisc[current_region] = False
            irdisc &= (last_entry > last_completion[:,current_region])
            
            for initial_region in init_regions[irdisc]:
                
                if track_transitions:
                    record_transition_entry(time_index,initial_region,current_region,weight)
                if track_eds:
                    ed = time_index - last_exit[initial_region]
                    record_ed_entry(time_index,initial_region,current_region,ed,weight)
                
                if track_fpts and last_completion[current_region,initial_region] > 0:
                        fpt = time_index - last_completion[current_region,initial_region]
                        record_fpt_entry(time_index,initial_region,current_region,fpt,weight)
                
                last_completion[initial_region,current_region] = time_index
                n_completions[initial_region,current_region] += 1
                         
            last_exit[last_region] = time_index
            last_entry[current_region] = time_index
            last_crossing[last_region,current_region] = time_index
                    
        self.time_index = tracking_table[-1,TIMEPOINT]
        self.last_region = tracking_table[-1,CURRENT_REGION]

            
def report_counts_averages(transacc, output_file):
    print('Number of crossings:',file=output_file)
    print(transacc.n_crossings,file=output_file)
    print('Number of completed transitions:',file=output_file)
    print(transacc.n_completions,file=output_file)
    
    if transacc.track_lifetimes:
        print('Average lifetimes:',file=output_file)
        print(transacc.lt_acc.average(),file=output_file)
        print('Standard deviation of dwell times:',file=output_file)
        print(transacc.lt_acc.std(),file=output_file)
    
    if transacc.track_eds:
        print('Average event durations:',file=output_file)
        print(transacc.ed_acc.average(),file=output_file)
        print('Standard deviation of event durations:',file=output_file)
        print(transacc.ed_acc.std(),file=output_file)

    if transacc.track_fpts:    
        print('Average first passage times:',file=output_file)
        print(transacc.fpt_acc.average(),file=output_file)
        print('Standard deviation of first passage times:',file=output_file)
        print(transacc.fpt_acc.std(),file=output_file)
    
    
def run_transanl_bf(args, transacc, pcoords):
    maxidx = args.stop or len(pcoords)
    n_chunks = int(math.ceil(maxidx/args.chunksize))
    for ichunk in xrange(0, n_chunks):
        lbi = ichunk*args.chunksize
        ubi = min((ichunk+1)*args.chunksize, maxidx)
        pc_chunk = pcoords[lbi:ubi]
        if pc_chunk.ndim == 1:
            pc_chunk = numpy.expand_dims(pc_chunk,-1)

        if ichunk > 0:    
            args.output_file.write('{current_time:d}/{total_time:d}\n'.format(current_time=lbi,
                                                                              total_time=maxidx))
            report_counts_averages(transacc, args.output_file)
            args.output_file.write('\n')
            args.output_file.flush()

        if ichunk == 0:
            transacc.accumulate_transitions(pc_chunk, continuation=False)
        else:
            transacc.accumulate_transitions(pc_chunk, continuation=True)
            
parser = wemd.rc.common_arg_parser('w_ttimes', description='''
Perform lifetime and transition analysis on WEMD or brute force data.''')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
#parser.add_argument('-I', '--save-converted-input', )
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
    print('loading bin data from WEMD system', file=args.output_file)
    runtime_config = wemd.rc.read_config(args.run_config_file)
    runtime_config.update_from_object(args)
    sim_manager = wemd.rc.load_sim_manager(runtime_config)
    sim_manager.load_data_manager()
    sim_manager.data_manager.open_backing()
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
    print('reading data from WEMD simulation as specified in', args.run_config_file, file=args.output_file)
    
args.output_file.write('Final results:\n')
report_counts_averages(transacc, args.output_file)
    