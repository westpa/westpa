from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math, warnings, re
import numpy, h5py
import wemd, wemdtools

import logging
log = logging.getLogger('w_ttimes')

from wemdtools.aframe import WEMDAnalysisTool,BinningMixin,WEMDDataReaderMixin,IterRangeMixin,MCBSMixin,TransitionAnalysisMixin

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

class WTTimes(MCBSMixin,TransitionAnalysisMixin,BinningMixin,IterRangeMixin,WEMDDataReaderMixin,WEMDAnalysisTool):
    def __init__(self):
        super(WTTimes,self).__init__()
        
        self.dt = None
        
        self.ed_stats_filename = None
        self.flux_stats_filename = None
        self.rate_stats_filename = None
        self.suppress_headers = None
        self.print_bin_labels = None

        self.ttimes_group = None        

        self.durations = None
        self.fluxes = None
        self.rates = None
        self.initial_bins = None
        self.final_bins = None

    def gen_stats(self):
        wemd.rc.pstatus('Analyzing transition statistics...')

        dt = self.dt
        n_sets = self.mcbs_nsets
        n_iters = self.last_iter - self.first_iter + 1
        lbi, ubi = self.calc_ci_bound_indices()
        pcoord_len = self.get_pcoord_len(self.first_iter)        
        total_time = n_iters * (pcoord_len - 1) * dt 
                
        transdat_ds = self.trans_h5group['transitions']
        transdat_ibin = transdat_ds['initial_bin']
        transdat_nblock = transdat_ds['block']
        transdat_in_range = (transdat_nblock >= self.first_iter) & (transdat_nblock <= self.last_iter) 
        
        durations = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        fluxes    = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        rates     = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        
        syn_avg_durations = numpy.empty((n_sets,), numpy.float64)
        syn_avg_fluxes    = numpy.empty((n_sets,), numpy.float64)
        syn_avg_rates     = numpy.empty((n_sets,), numpy.float64)
        
        w_n_bins = len(str(self.n_bins))
        w_n_sets = len(str(n_sets))
        
        for ibin in self.initial_bins:
            trans_ibin = transdat_ds[(transdat_ibin == ibin) & transdat_in_range]
            for fbin in self.final_bins:
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
                    wemd.rc.pstatus('\r  {:{w_n_bins}d}->{:<{w_n_bins}d} set {:{w_n_sets}d}/{:<{w_n_sets}d}, set size {:<20d}'\
                            .format(ibin,fbin,iset+1,n_sets,dlen, w_n_bins=w_n_bins, w_n_sets=w_n_sets), end='')
                    wemd.rc.pflush()
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
            wemd.rc.pstatus()
            del trans_ibin
            
        for (dsname, data) in (('duration', durations), ('flux', fluxes), ('rate', rates)):
            try:
                del self.ttimes_group[dsname]
            except KeyError:
                pass
            
            ds = self.ttimes_group.create_dataset(dsname, data=data)
            attrs = ds.attrs
            attrs['dt'] = dt
            attrs['total_time'] = total_time
            attrs['ci_alpha'] = self.mcbs_alpha
            attrs['ci_n_sets'] = self.mcbs_nsets
            
            self.record_data_iter_range(ds)
            self.record_data_binhash(ds)

        attrs = self.ttimes_group.attrs
        attrs = ds.attrs
        attrs['dt'] = dt
        attrs['total_time'] = total_time
        attrs['ci_alpha'] = self.mcbs_alpha
        attrs['ci_n_sets'] = self.mcbs_nsets
        
            
        self.durations = durations
        self.fluxes = fluxes
        self.rates = rates
        
        self.record_data_iter_range(self.ttimes_group)
        self.record_data_binhash(self.ttimes_group)
            
    def summarize_stats(self):
        for (array, dsname, argname, title) in ((self.durations, 'duration', 'ed_stats', 'event duration'),
                                                (self.fluxes, 'flux', 'flux_stats', 'flux'),
                                                (self.rates, 'rate', 'rate_stats', 'rate')):
            if array is None:
                array = self.ttimes_group[dsname]
            self.summarize_ci(getattr(args,argname), array, title, self.mcbs_display_confidence,
                              headers=(not self.suppress_headers), labels=self.print_bin_labels)
        
            
    def summarize_ci(self, filename, array, title, confidence, headers, labels):
        if not filename: return
        
        format_2d = '{ibin:{mw}d}    {fbin:{mw}d}    {0:20.15g}    {1:20.15g}    {2:20.15g}    {3:20.15g}    {4:20.15g}    {5:20.15g}\n'
        max_ibin_width = len(str(self.n_bins-1))
         
        outfile = open(filename, 'wt')
        if headers:
            outfile.write('''\
# {title:} statistics
# confidence interval = {confidence}%
# ----
# column 0: initial bin index
# column 1: final bin index
# column 2: lower bound of confidence interval
# column 3: upper bound of confidence interval
# column 4: width of confidence interval
# column 5: relative width of confidence interval [abs(width/average)]
# column 6: symmetrized error [max(upper-average, average-lower)]
# ----
'''.format(title=title, confidence=confidence))
            if labels:
                self.write_bin_labels(outfile)
                outfile.write('----\n')

        for ibin in self.initial_bins:
            for fbin in self.final_bins:
                mean = array[ibin,fbin]['expectation']
                lb = array[ibin,fbin]['ci_lower']
                ub = array[ibin,fbin]['ci_upper']
                ciwidth = ub - lb
                relciwidth = abs(ciwidth/mean)
                symmerr = max(mean-lb,ub-mean)
                
                outfile.write(format_2d.format(*map(float,(mean,lb,ub,ciwidth,relciwidth,symmerr)),
                                               ibin=ibin,fbin=fbin,mw=max_ibin_width))
     
    def parse_simple_int_range(self, range_string):
        try:
            entries = []
            fields = re.split('\s*,\s*', range_string)
            for field in fields:
                if '-' in field:
                    lb, ub = map(int,re.split('\s*-\s*', field))
                    entries.extend(range(lb,ub+1))
                else:
                    entries.append(int(field))
        except (ValueError,TypeError):
            raise ValueError('invalid range string {!r}'.format(range_string))
        else:
            return entries

        

###################################################################################################################################

wtt = WTTimes()
parser = argparse.ArgumentParser('w_ttimes', description='''\
Trace the WEMD trajectory tree and report on transition kinetics.
''')
wemd.rc.add_args(parser)
wtt.add_args(parser)

agroup = parser.add_argument_group('kinetics analysis options')
agroup.add_argument('--dt', dest='dt', type=float, default=1.0,
                    help='Assume input data has a time spacing of DT (default: %(default)s).')
agroup.add_argument('-i', '--initial-bins', dest='ibins_string', metavar='I',
                     help='''Only calculate statistics for transitions starting in bin I.  This may be specified as a
                     comma-separated list of integers or ranges, as in "0,2-4,5,9"''')
agroup.add_argument('-j', '--final-bins', dest='fbins_string', metavar='J',
                    help='''Only calculate statistics for transitions ending in bin J.  This may be specified as a
                     comma-separated list of integers or ranges, as in "0,2-4,5,9"''')

output_options = parser.add_argument_group('kinetics analysis output options')        
output_options.add_argument('--edstats', dest='ed_stats', default='edstats.txt',
                            help='Store event duration statistics in ED_STATS (default: edstats.txt)')
output_options.add_argument('--fluxstats', dest='flux_stats', default='fluxstats.txt',
                            help='Store flux statistics in FLUX_STATS (default: fluxstats.txt)')
output_options.add_argument('--ratestats', dest='rate_stats', default='ratestats.txt',
                            help='Store rate statistics in RATE_STATS (default: ratestats.txt)')
output_options.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                            help='Do not include headers in text output files (default: include headers)')
output_options.add_argument('--binlabels', dest='print_bin_labels', action='store_true',
                            help='Print bin labels in output files, if available (default: do not print bin labels)')


args = parser.parse_args()
wemd.rc.process_args(args, config_required = False)
wtt.process_args(args)

wtt.dt = args.dt
wtt.ed_stats_filename = args.ed_stats
wtt.flux_stats_filename = args.flux_stats
wtt.rate_stats_filename = args.rate_stats
wtt.suppress_headers = args.suppress_headers
wtt.print_bin_labels = args.print_bin_labels

if args.ibins_string:
    wtt.initial_bins = wtt.parse_simple_int_range(args.ibins_string)
    wemd.rc.pstatus('Will report statistics for transitions beginning in the following bins:\n{!s}'.format(wtt.initial_bins))
else:
    wemd.rc.pstatus('Will report statistics for transitions beginning in any bin.')
    wtt.initial_bins = range(wtt.n_bins)
    
if args.fbins_string:
    wtt.final_bins = wtt.parse_simple_int_range(args.fbins_string)
    wemd.rc.pstatus('Will report statistics for transitions ending in the following bins:\n{!s}'.format(wtt.final_bins))
else:
    wemd.rc.pstatus('Will report statistics for transitions ending in any bin.')
    wtt.final_bins = range(wtt.n_bins) 

wtt.check_iter_range()
wtt.open_analysis_backing()
wtt.ttimes_group = wtt.require_analysis_group('w_ttimes', replace=False)
wtt.require_bin_assignments()
wtt.require_transitions()
wtt.gen_stats()
wtt.summarize_stats()

