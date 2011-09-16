from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math, warnings, re
import numpy, h5py
import wemd, wemdtools
from wemdtools.miscfn import parse_int_list
from wemdtools.files import load_npy_or_text

import logging
log = logging.getLogger('w_ttimes_bf')

from wemdtools.aframe import WEMDAnalysisTool,BinningMixin,MCBSMixin,TransitionAnalysisMixin,KineticsAnalysisMixin
from wemdtools.aframe.transitions import TransitionEventAccumulator

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

class WTTimesBF(KineticsAnalysisMixin,MCBSMixin,TransitionAnalysisMixin,BinningMixin,WEMDAnalysisTool):
    def __init__(self):
        super(WTTimesBF,self).__init__()
                
        self.ed_stats_filename = None
        self.fpt_stats_filename = None
        self.suppress_headers = None
        self.print_bin_labels = None
        
        self.chunksize = None
        self.usecols = None
        self.input_files = None

        self.ttimes_group = None        

        self.durations = None
        self.fpts = None

    def check_bin_data(self):
        '''Check to see that existing binning data corresponds to the same bin topology and iteration range as requested'''
        
        self.require_binning_group()
        
        if self.discard_bin_assignments:
            #wemd.rc.pstatus('Discarding existing binning data.')
            self.delete_binning_group()
        elif 'bin_assignments' in self.binning_h5group:
            if not self.check_data_binhash(self.binning_h5group):
                wemd.rc.pstatus('Bin definitions have changed; deleting existing binning data.')
                self.delete_binning_group()
        else:
            wemd.rc.pstatus('Using existing bin assignments.')
        self.require_binning_group()
        
    def check_transitions(self):
        self.require_transitions_group()
        
        if self.discard_transition_data:
            wemd.rc.pstatus('Discarding existing transition data.')
            self.delete_transitions_group()
        elif not self.check_data_binhash(self.trans_h5group):
            wemd.rc.pstatus('Bin definitions have changed; deleting existing transition data.')
            self.delete_transitions_group()
        else:
            wemd.rc.pstatus('Using existing transition data.')
                
        self.require_transitions_group()

    def assign_to_bins(self):
        wemd.rc.pstatus('Assigning to bins...')
        
        assignments_ds = self.binning_h5group.create_dataset('bin_assignments', dtype=numpy.min_scalar_type(self.n_bins),
                                                             shape=(1,), maxshape=(None,), chunks=(65536,),
                                                             compression='gzip',)
        
        trajlen = 0
        for (ifile, filename) in enumerate(self.input_files):
            if ifile > 0:
                wemd.rc.pstatus()
            input_data = load_npy_or_text(filename)
            nrows = len(input_data)
            maxlen_nrows = len(str(nrows))
            pct_prec = max(0,int(math.ceil(-math.log10(self.chunksize/nrows)))-2)
            
            for istart in xrange(0,nrows,self.chunksize):
                if self.usecols:
                    pcoords = input_data[istart:istart+self.chunksize,self.usecols]
                else:
                    pcoords = input_data[istart:istart+self.chunksize]
                    
                assignments = self.region_set.map_to_all_indices(pcoords)
                assignments_ds.resize((trajlen+len(assignments),))
                assignments_ds[trajlen:trajlen+len(assignments)] = assignments
                trajlen += len(assignments)
                wemd.rc.pstatus('\r  {:s}: {:{mlnr}d}/{:<{mlnr}d} = {:.{pp}f}%'.format(filename, istart+len(assignments), nrows,
                                                                                      (istart+len(assignments))/nrows*100,
                                                                                      mlnr=maxlen_nrows,pp=pct_prec), end='')
                wemd.rc.pflush()
                del pcoords
            del input_data
        wemd.rc.pstatus()
        for h5object in (self.binning_h5group, assignments_ds):
            self.record_data_binhash(h5object)
                
        wemd.rc.pstatus()
            
    def require_bin_assignments(self):
        self.check_bin_data()
                
        # The group will be empty if the user requested the data to go away, or if the data is not conformant
        # with the current bins and iteration range, so the following lets us know if we need
        # to recalculate 
        if not 'bin_assignments' in self.binning_h5group:
            self.assign_to_bins()

    def find_transitions(self):
        wemd.rc.pstatus('Finding transitions...')
        self.require_binning_group()
        output_group = self.require_analysis_group('transitions')
            
        self.accumulator = TransitionEventAccumulator(self.n_bins, output_group, calc_fpts = False)
        assignments_ds = self.binning_h5group['bin_assignments']
        
        nrows = assignments_ds.len()
        chunksize = int(self.chunksize / 10)
        maxwidth_nrows = len(str(nrows))
        for istart in xrange(0, nrows, chunksize):
            assignments = assignments_ds[istart:min(istart+chunksize,nrows)]
            weights = numpy.ones((len(assignments),))
            binpops = numpy.ones((len(assignments),self.n_bins))
            
            if istart == 0:
                self.accumulator.start_accumulation(assignments, weights, binpops)
            else:
                self.accumulator.continue_accumulation(assignments, weights, binpops)
                
            wemd.rc.pstatus('\r  {:{mwnr}d}/{:<{mwnr}d} ({:.2f}%)'.format(istart+len(assignments), nrows, 
                                                                       (istart+len(assignments))/nrows*100, mwnr=maxwidth_nrows), 
                            end='')
            wemd.rc.pflush()
            
            self.accumulator.flush_transition_data()
            del assignments, weights, binpops
        wemd.rc.pstatus()                 
                
        try:
            del output_group['n_trans']
        except KeyError:
            pass
        output_group['n_trans'] = self.accumulator.n_trans
        
        for h5object in (output_group, output_group['n_trans'], output_group['transitions']):
            self.record_data_binhash(h5object)
        
        self.accumulator.clear()
    

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
            ds.attrs.update({'dt': dt, 'total_time': total_time, 'ci_alpha': self.mcbs_alpha, 'ci_n_sets': self.mcbs_nsets})
            
            self.record_data_iter_range(ds)
            self.record_data_binhash(ds)

        self.ttimes_group.attrs.update({'dt': dt, 'total_time': total_time, 
                                        'ci_alpha': self.mcbs_alpha, 'ci_n_sets': self.mcbs_nsets})
            
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
     

        

###################################################################################################################################

wtt = WTTimesBF()
parser = argparse.ArgumentParser('w_ttimes_bf', description='''\
Calculate transition kinetics information for brute force simulations.
''')
wemd.rc.add_args(parser)
wtt.add_args(parser)
parser.set_defaults(anal_h5name='bf_analysis.h5')

output_options = parser.add_argument_group('kinetics analysis output options')        
output_options.add_argument('--edstats', dest='ed_stats', default='edstats_bf.txt',
                            help='Store event duration statistics in ED_STATS (default: edstats_bf.txt)')
output_options.add_argument('--fptstats', dest='fpt_stats', default='fptstats_bf.txt',
                            help='Store event duration statistics in ED_STATS (default: fptstats_bf.txt)')
output_options.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                            help='Do not include headers in text output files (default: include headers)')
output_options.add_argument('--binlabels', dest='print_bin_labels', action='store_true',
                            help='Print bin labels in output files, if available (default: do not print bin labels)')

input_options = parser.add_argument_group('input options')
input_options.add_argument('datafiles', nargs='+', metavar='DATAFILE',
                           help='''Trajectory file(s) to analyze, either text or Numpy (.npy or .npz) format.  Uncompressed numpy
                           files will be memory-mapped, allowing analysis of data larger than available RAM (though not
                           larger than the available address space).''')
input_options.add_argument('--usecols', dest='usecols', metavar='COLUMNS', type=parse_int_list,
                           help='''Use only the given COLUMNS from the input file(s), e.g. "0", "0,1", "0:5,7,9:10".''')
input_options.add_argument('--chunksize', dest='chunksize', type=long, default=2097152,
                           help='''Process input data in blocks of size CHUNKSIZE.  This will only reduce memory requirements
                           when using Numpy (.npy) format input. (Default: %(default)d.)''')

args = parser.parse_args()
wemd.rc.process_args(args, config_required = False)
wtt.process_args(args)

if args.usecols:
    wemd.rc.pstatus('Using only the following columns from the input file: {!s}'.format(args.usecols))
    wtt.usecols = args.usecols
else:
    wtt.usecols = None
    
wtt.ed_stats_filename = args.ed_stats
wtt.fpt_stats_filename = args.fpt_stats
wtt.suppress_headers = args.suppress_headers
wtt.print_bin_labels = args.print_bin_labels
wtt.input_files = args.datafiles
wtt.chunksize = args.chunksize

wtt.open_analysis_backing()
wtt.ttimes_group = wtt.require_analysis_group('w_ttimes_bf', replace=False)
wtt.check_bin_data()
wtt.check_transitions()
wtt.require_bin_assignments()
wtt.require_transitions()
#wtt.gen_stats()
#wtt.summarize_stats()

