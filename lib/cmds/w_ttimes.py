from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math, warnings
import numpy, h5py
import wemd, wemdtools
#from wemdtools.transitions.transacc import TransitionEventAccumulator
#from wemdtools.trajectories.trajtree import TrajTree

import logging
log = logging.getLogger('w_ttimes')

from wemdtools.aframe import WEMDAnalysisTool,BinningMixin,DataReaderMixin,IterRangeMixin,MCBSMixin,TransitionAnalysisMixin

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

class WTTimes(MCBSMixin,TransitionAnalysisMixin,BinningMixin,IterRangeMixin,DataReaderMixin,WEMDAnalysisTool):
    def __init__(self):
        super(WTTimes,self).__init__()
        
        self.dt = None
        self.suppress_headers = None
        self.wtt_group = None
        
        self.trajtree = None
        self.accumulator = None
        
        self.assignments = None
        self.populations = None
        self.durations = None
        self.fluxes = None
        self.rates = None
        
        self.ttimes_group = None
        

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

wtt = WTTimes()
parser = argparse.ArgumentParser('w_ttimes', description='''\
Trace the WEMD trajectory tree and report on transition kinetics.
''')
wemd.rc.add_args(parser)
wtt.add_args(parser)

args = parser.parse_args()
wemd.rc.process_args(args, config_required = False)
wtt.process_args(args)

wtt.check_iter_range()
wtt.open_analysis_backing()
wtt.ttimes_group = wtt.require_analysis_group('w_ttimes', replace=False)
wtt.require_bin_assignments()
wtt.require_transitions()



#try:
#    h5fn, h5gn = args.h5.rsplit(':',1)
#except ValueError:
#    # No colon
#    h5fn = args.h5
#    h5gn = 'w_ttimes'
#    
#h5file = h5py.File(h5fn)
#h5group = h5file.require_group(h5gn)
#    
#sim_manager = wemdtools.files.sim_manager_from_args(args, 'datafile', status_stream=sys.stdout)
#data_manager = wemdtools.data_manager.CachingWEMDDataReader(sim_manager.load_data_manager(), cache_pcoords=False)
#region_set = wemdtools.bins.get_region_set_from_args(args, status_stream=sys.stdout)
#
#first_iter = args.first_iter
#last_iter = args.last_iter or data_manager.current_iteration - 1
#last_iter = min(last_iter, data_manager.current_iteration - 1)
#
#helper = TTimesHelper(data_manager, region_set, first_iter, last_iter, h5file, h5group, args.quiet_mode)
#n_iters = helper.n_iters
#n_bins = helper.n_bins
#
#if n_iters == 0:
#    sys.stdout.write('No data to analyze; exiting.  If this seems wrong, check --start and --stop.\n')
#    sys.exit(0)
#    
#if args.limit_calc_to in (None,'bins'):
#    helper.assign_to_bins()
#
#if args.limit_calc_to in (None, 'transitions'):
#    helper.scan_trajectories()
#
#if args.limit_calc_to in (None, 'statistics'):
#    alpha = 1.0 - args.confidence
#    n_sets = args.bssize or wemdtools.stats.mcbs.get_bssize(alpha)    
#    helper.gen_stats(n_sets, alpha, args.dt)    
#    
#if args.limit_calc_to in (None, 'summary'):
#    helper.summarize_stats(args)
        
