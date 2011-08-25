from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math, warnings
import numpy, scipy, h5py
import scipy.stats
import wemd, wemdtools
from wemdtools.transitions.transacc import TransitionEventAccumulator

import logging
log = logging.getLogger('w_ttimes')

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

def pstatus(*fnargs, **fnkwargs):
    global args
    if not args.quiet_mode:
        outfile = fnkwargs.get('file', sys.stdout)
        print(*fnargs,**fnkwargs)
        try:
            outfile.flush()
        except IOError:
            pass    


class TTimesHelper:
    def __init__(self, data_manager, region_set, first_iter, last_iter, output_h5file, output_h5group, quiet_mode = False):
        self.sim_manager = sim_manager
        self.args = args        
        self.output_h5file = output_h5file
        self.output_h5group = output_h5group
        
        self.region_set = region_set
        self.n_bins = len(region_set.get_all_bins())
        
        self.data_manager = data_manager
        self.trajtree = wemdtools.trajectories.trajtree.TrajTree(self.data_manager, include_pcoords = False)
        self.accumulator = wemdtools.transitions.transacc.TransitionEventAccumulator(self.n_bins, self.output_h5group)
        
        self.first_iter = first_iter
        self.last_iter = last_iter
        self.n_iters = last_iter - first_iter + 1
        
        n_particles = data_manager.h5file['/summary']['n_particles']
        self.n_total_segs = n_particles[self.first_iter-1:self.last_iter].sum() 
        self.n_segs_max   = n_particles[self.first_iter-1:self.last_iter].max()
        self.pcoord_len = data_manager.get_iter_group(self.first_iter)['pcoord'].shape[1]

        self.assignments = None
        self.populations = None
        self.durations = None
        self.fluxes = None
        self.rates = None
        
        self.quiet_mode = quiet_mode
        
    
    def assign_to_bins(self):
        pstatus('Assigning to bins...')
        for key in ('bin_assignments', 'bin_populations'):
            try:
                del self.output_h5group[key]
            except KeyError:
                pass
        
        assignments = numpy.zeros((n_iters, self.n_segs_max, self.pcoord_len), dtype=numpy.min_scalar_type(n_bins))
        populations = numpy.zeros((n_iters, self.pcoord_len, n_bins), dtype=numpy.float64)
        
        for (iiter, n_iter) in enumerate(xrange(first_iter, last_iter+1)):
            pstatus('\r  Iteration {:d}'.format(n_iter), end='')
            seg_index = data_manager.get_seg_index(n_iter) 
            pcoords = data_manager.get_pcoord_array(n_iter) # pcoords is [seg_id, time, dimension]
            weights = seg_index['weight'] # weight is [seg_id]
            
            for seg_id in xrange(len(seg_index)):
                assignments[iiter,seg_id,:] = region_set.map_to_all_indices(pcoords[seg_id,:,:]) # results in time-resolved indices
            
            for it in xrange(self.pcoord_len):
                populations[iiter, it, :] = numpy.bincount(assignments[iiter,:len(seg_index),it], weights, minlength=n_bins)
                        
            del pcoords, weights, seg_index
        
        # Compress bin assignments, as they are infrequently used and contain lots of zero bits
        self.output_h5group.create_dataset('bin_assignments', data=assignments, compression='gzip')
        self.output_h5group['bin_assignments'].attrs['first_iter'] = self.first_iter
        
        # Might as well compress populations, too, since we just read the whole dang thing before
        # we scan trajectories anyway
        self.output_h5group.create_dataset('bin_populations', data=populations, compression='gzip')
        self.output_h5group['bin_populations'].attrs['first_iter'] = self.first_iter
        
        self.assignments = assignments
        self.populations = populations
        pstatus()
        
    def scan_trajectories(self):
        pstatus('Finding transitions...')
        if self.assignments is None:
            self.assignments = self.output_h5group['bin_assignments'][...]
        if self.populations is None:
            self.populations = self.output_h5group['bin_populations'][...]
            
        self.n_segs_visited = 0
        self.trajtree.trace_trajectories(self.first_iter, self.last_iter, callable=self._segment_callback)
        self.accumulator.flush_transition_data()
        self.output_h5group['n_trans'] = self.accumulator.n_trans
        self.accumulator.clear()
        self.data_manager.clear_cache()
        pstatus()
        
        del self.assignments, self.populations
        self.assignments = self.populations = None
        
    def _segment_callback(self, segment, children, history):
        iiter = segment.n_iter - self.first_iter
        seg_id = segment.seg_id
        weights = numpy.empty((self.pcoord_len,), numpy.float64)
        weights[:] = segment.weight
        bin_pops = self.populations[iiter, :, :]
        
        if len(history) == 0:
            # New trajectory
            self.accumulator.start_accumulation(self.assignments[iiter, seg_id, :], weights, bin_pops, block=segment.n_iter)
        else:
            # Continuing trajectory
            self.accumulator.continue_accumulation(self.assignments[iiter, seg_id, :], weights, bin_pops, block=segment.n_iter)
            
        self.n_segs_visited += 1
        
        if not self.quiet_mode and (self.n_segs_visited % 1000 == 0 or self.n_segs_visited == self.n_total_segs):
            pct_visited = self.n_segs_visited / self.n_total_segs * 100
            pstatus('\r  {:d} of {:d} segments ({:.1f}%) analyzed'.format(long(self.n_segs_visited), long(self.n_total_segs), 
                                                                          float(pct_visited)), end='')

    def gen_stats(self, n_sets, alpha, dt):
        pstatus('Analyzing transition statistics...')

        lbi = long(math.floor(n_sets*alpha/2))
        ubi = long(math.ceil(n_sets*(1-alpha/2)))
        
        total_time = self.n_iters * (self.pcoord_len - 1) * dt 
        
        pstatus('  Using bootstrap sampling of {:d} sets'.format(n_sets))
        pstatus('  {:g}% confidence interval (alpha={:g}) requested [bounding indices ({:d},{:d})]'\
                .format(args.confidence*100, alpha, lbi, ubi))
        
        transdat_ds = self.output_h5group['transitions']
        transdat_ibin = transdat_ds['initial_bin']
        
        durations = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        fluxes    = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        rates     = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        
        syn_avg_durations = numpy.empty((n_sets,), numpy.float64)
        syn_avg_fluxes    = numpy.empty((n_sets,), numpy.float64)
        syn_avg_rates     = numpy.empty((n_sets,), numpy.float64)
        
        w_n_bins = len(str(n_bins))
        w_n_sets = len(str(n_sets))
        
        for ibin in xrange(self.n_bins):
            trans_ibin = transdat_ds[transdat_ibin == ibin]
            for fbin in xrange(self.n_bins):
                trans_ifbins = trans_ibin[trans_ibin['final_bin'] == fbin]
                dlen = len(trans_ifbins)
                
                if not dlen: continue
                
                trans_weights = trans_ifbins['final_weight']
                trans_durations = trans_ifbins['duration']
                trans_ibinprobs = trans_ifbins['initial_bin_pop']
                                    
                durations[ibin,fbin]['expectation'] = numpy.average(trans_durations, weights=trans_weights) * dt
                avg_flux = trans_weights.sum() / total_time
                fluxes[ibin,fbin]['expectation']    = avg_flux
                rates[ibin,fbin]['expectation']     = avg_flux / trans_ibinprobs.mean()
                
                for iset in xrange(n_sets):
                    pstatus('\r  {:{w_n_bins}d}->{:<{w_n_bins}d} set {:{w_n_sets}d}/{:<{w_n_sets}d}, set size {:d}'\
                            .format(ibin,fbin,iset+1,n_sets,dlen, w_n_bins=w_n_bins, w_n_sets=w_n_sets), end='')
                    indices = numpy.random.randint(dlen, size=(dlen,))
                    syn_weights   = trans_weights[indices]
                    syn_durations = trans_durations[indices]
                    syn_ibinprobs = trans_ibinprobs[indices]
                    
                    syn_avg_durations[iset] = numpy.average(syn_durations, weights=syn_weights) * dt
                    syn_avg_fluxes[iset] = syn_weights.sum() / total_time
                    syn_avg_rates[iset] = syn_avg_fluxes[iset] / syn_ibinprobs.mean()
                    
                    del indices, syn_weights, syn_durations, syn_ibinprobs
                    
                syn_avg_durations.sort()
                syn_avg_fluxes.sort()
                syn_avg_rates.sort()
                
                durations[ibin,fbin]['ci_lower'] = syn_avg_durations[lbi]
                durations[ibin,fbin]['ci_upper'] = syn_avg_durations[ubi]
                
                fluxes[ibin,fbin]['ci_lower'] = syn_avg_fluxes[lbi]
                fluxes[ibin,fbin]['ci_upper'] = syn_avg_fluxes[ubi]
                
                rates[ibin,fbin]['ci_lower'] = syn_avg_rates[lbi]
                rates[ibin,fbin]['ci_upper'] = syn_avg_rates[ubi]
                    
                del trans_weights, trans_durations, trans_ibinprobs, trans_ifbins
            pstatus()
            del trans_ibin
            
        for (dsname, data) in (('duration', durations), ('flux', fluxes), ('rate', rates)):
            try:
                del self.output_h5group[dsname]
            except KeyError:
                pass
            
            self.output_h5group[dsname] = data
            self.output_h5group[dsname].attrs['dt'] = dt
            self.output_h5group[dsname].attrs['total_time'] = total_time
            self.output_h5group[dsname].attrs['alpha'] = alpha
            
        self.durations = durations
        self.fluxes = fluxes
        self.rates = rates
            
    def summarize_stats(self, args):
        for (array, dsname, argname, title) in ((self.durations, 'duration', 'ed_stats', 'event duration'),
                                                (self.fluxes, 'flux', 'flux_stats', 'flux'),
                                                (self.rates, 'rate', 'rate_stats', 'rate')):
            if array is None:
                array = self.output_h5group[dsname]
                self.summarize_ci(getattr(args,argname), array, title, args.confidence,
                                  headers=(not args.suppress_headers), labels=args.print_labels)
        
            
    def summarize_ci(self, filename, array, title, confidence, headers, labels):
        if not filename: return
        
        format_2d = '{ibin:{mw}d}    {fbin:{mw}d}    {0:20.15g}    {1:20.15g}    {2:20.15g}    {3:20.15g}    {4:20.15g}    {5:20.15g}\n'
        max_ibin_width = len(str(self.n_bins-1))
         
        outfile = open(filename, 'wt')
        if headers:
            outfile.write('''\
# {title:} statistics
# confidence interval = {confidence:g}%
# ----
# column 0: initial bin index
# column 1: final bin index
# column 2: lower bound of confidence interval
# column 3: upper bound of confidence interval
# column 4: width of confidence interval
# column 5: relative width of confidence interval [abs(width/average)]
# column 6: symmetrized error [max(upper-average, average-lower)]
# ----
'''.format(title=title, confidence=confidence*100))
            if labels:
                wemdtools.bins.print_labels(self.region_set, outfile)
                outfile.write('----\n')

        for ibin in xrange(n_bins):
            for fbin in xrange(n_bins):
                mean = array[ibin,fbin]['expectation']
                lb = array[ibin,fbin]['ci_lower']
                ub = array[ibin,fbin]['ci_upper']
                ciwidth = ub - lb
                relciwidth = abs(ciwidth/mean)
                symmerr = max(mean-lb,ub-mean)
                
                outfile.write(format_2d.format(*map(float,(mean,lb,ub,ciwidth,relciwidth,symmerr)),
                                               ibin=ibin,fbin=fbin,mw=max_ibin_width))
     
    

        

###################################################################################################################################        
parser = wemd.rc.common_arg_parser('w_ttimes', 'Kinetics analysis for WEMD data')
parser.add_argument('-q', '--quiet', dest='quiet_mode', action='store_true',
                    help='Do not emit status messages.')
input_options = parser.add_argument_group('input')
input_options.add_argument('datafile', nargs='?', metavar='DATAFILE',
                           help='Input HDF5 file. (Default: use HDF5 specified in wemd.cfg')
input_options.add_argument('-b', '--begin', '--start', '--first', dest='first_iter', type=int, default=1,
                           help='Begin processing at iteration START (default: 1)')
input_options.add_argument('-e', '--end', '--stop', '--last', dest='last_iter', type=int,
                           help='Stop processing after iteration STOP (default: last completed iteration).')


calc_options = parser.add_argument_group('calculation and intermediate storage')
calc_options.add_argument('-H', '--h5', metavar='H5FILE:H5GROUP', default='analysis.h5:/w_ttimes',
                          help='Store intermediate values and results in H5FILE, in group H5GROUP. '
                          +'(Default: analysis.h5:/w_ttimes)')
calc_options.add_argument('--dt', type=float, default=1.0,
                          help='Assume input data has a time spacing of DT (default: 1.0)')

only_options = calc_options.add_mutually_exclusive_group(required=False)
only_options.add_argument('--assign-only', dest='limit_calc_to', action='store_const', const='bins',
                          help='Perform bin assignments only.')
only_options.add_argument('--transitions-only', dest='limit_calc_to', action='store_const', const='transitions',
                          help='Perform transition analysis only.')
only_options.add_argument('--statistics-only', dest='limit_calc_to', action='store_const', const='statistics',
                          help='Perform statistical analysis on transitions and kinetics only.')
only_options.add_argument('--summarize-only', dest='limit_calc_to', action='store_const', const='summary',
                          help='Only write statistics from previous analysis.')

output_options = parser.add_argument_group('output options')
output_options.add_argument('-E', '--edstats', dest='ed_stats', default='edstats.txt',
                            help='Store event duration statistics in ED_STATS (default: edstats.txt)')
output_options.add_argument('-f', '--fluxstats', dest='flux_stats', default='fluxstats.txt',
                            help='Store flux statistics in FLUX_STATS (default: fluxstats.txt)')
output_options.add_argument('-R', '--ratestats', dest='rate_stats', default='ratestats.txt',
                            help='Store rate statistics in RATE_STATS (default: ratestats.txt)')
output_options.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                            help='Do not include headers in text output files (default: include headers)')
output_options.add_argument('-l', '--labels', dest='print_labels', action='store_true',
                            help='Print bin labels in output files, if available (default: do not print bin labels)')

wemdtools.bins.add_region_set_options(parser)
wemdtools.stats.mcbs.add_mcbs_options(parser)

args = parser.parse_args()

try:
    h5fn, h5gn = args.h5.rsplit(':',1)
except ValueError:
    # No colon
    h5fn = args.h5
    h5gn = 'w_ttimes'
    
h5file = h5py.File(h5fn)
h5group = h5file.require_group(h5gn)
    
sim_manager = wemdtools.files.sim_manager_from_args(args, 'datafile', status_stream=sys.stdout)
data_manager = wemdtools.data_manager.CachingDataReader(sim_manager.load_data_manager(), cache_pcoords=False)
region_set = wemdtools.bins.get_region_set_from_args(args, status_stream=sys.stdout)

first_iter = args.first_iter
last_iter = args.last_iter or data_manager.current_iteration - 1
last_iter = min(last_iter, data_manager.current_iteration - 1)

helper = TTimesHelper(data_manager, region_set, first_iter, last_iter, h5file, h5group, args.quiet_mode)
n_iters = helper.n_iters
n_bins = helper.n_bins

if n_iters == 0:
    sys.stdout.write('No data to analyze; exiting.  If this seems wrong, check --start and --stop.\n')
    sys.exit(0)
    
if args.limit_calc_to in (None,'bins'):
    helper.assign_to_bins()

if args.limit_calc_to in (None, 'transitions'):
    helper.scan_trajectories()

if args.limit_calc_to in (None, 'statistics'):
    alpha = 1.0 - args.confidence
    n_sets = args.bssize or wemdtools.stats.mcbs.get_bssize(alpha)    
    helper.gen_stats(n_sets, alpha, args.dt)    
    
if args.limit_calc_to in (None, 'summary'):
    helper.summarize_stats(args)
        
