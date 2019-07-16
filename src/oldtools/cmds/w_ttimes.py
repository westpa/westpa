# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.


import argparse, math
import numpy
import westpa, oldtools

import logging
log = logging.getLogger('w_ttimes')

from oldtools.aframe import (WESTAnalysisTool,BinningMixin,WESTDataReaderMixin,IterRangeMixin,MCBSMixin,TransitionAnalysisMixin,
                              KineticsAnalysisMixin,CommonOutputMixin,BFDataManager,BFTransitionAnalysisMixin)

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

class WTTimesBase:
    def __init__(self):
        super(WTTimesBase,self).__init__()
        
        self.ed_stats_filename = None
        self.fpt_stats_filename = None
        self.flux_stats_filename = None
        self.rate_stats_filename = None
        self.suppress_headers = None
        self.print_bin_labels = None

        self.ttimes_group = None        

        self.durations = None
        self.fpts = None
        self.fluxes = None
        self.rates = None
        
        
    def add_args(self, parser, upcall = True):
        '''Add arguments to a parser common to all analyses of this type.'''
        if upcall:
            try:
                upfunc = super(WTTimesBase,self).add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)

        output_options = parser.add_argument_group('kinetics analysis output options')        
        output_options.add_argument('--edstats', dest='ed_stats', default='edstats.txt',
                                    help='Store event duration statistics in ED_STATS (default: edstats.txt)')
        if self.bf_mode:
            output_options.add_argument('--fptstats', dest='fpt_stats', default='fptstats.txt',
                                        help='Store first passage time statistics in FPT_STATS (default: fptstats.txt).')
        else:
            output_options.add_argument('--fptstats', dest='fpt_stats',
                                        help='Store first passage time statistics in FPT_STATS (default: do not store).')
        output_options.add_argument('--fluxstats', dest='flux_stats', default='fluxstats.txt',
                                    help='Store flux statistics in FLUX_STATS (default: fluxstats.txt)')
        output_options.add_argument('--ratestats', dest='rate_stats', default='ratestats.txt',
                                    help='Store rate statistics in RATE_STATS (default: ratestats.txt)')
        self.add_common_output_args(output_options)
        
    def process_args(self, args, upcall = True):
        
        self.ed_stats_filename = args.ed_stats
        self.fpt_stats_filename = args.fpt_stats
        self.flux_stats_filename = args.flux_stats
        self.rate_stats_filename = args.rate_stats
        self.process_common_output_args(args)
        self.calc_fpts = bool(args.fpt_stats)
        
        if upcall:
            try:
                upfunc = super(WTTimesBase,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)        
                
    def gen_stats(self):
        self.require_transitions_group()
        westpa.rc.pstatus('Analyzing transition statistics...')

        dt = self.dt
        n_sets = self.mcbs_nsets
        lbi, ubi = self.calc_ci_bound_indices()

        total_time = self.get_total_time()                
        transdat_ds = self.trans_h5group['transitions']
        transdat_ibin = transdat_ds['initial_bin']
        
        if not self.bf_mode:
            transdat_niter = transdat_ds['n_iter']
            transdat_in_range = (transdat_niter >= self.first_iter) & (transdat_niter <= self.last_iter) 
        
        durations = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        fpts      = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        fluxes    = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        rates     = numpy.zeros((self.n_bins,self.n_bins), ciinfo_dtype)
        
        syn_avg_durations = numpy.empty((n_sets,), numpy.float64)
        syn_avg_fpts      = numpy.empty((n_sets,), numpy.float64)
        syn_avg_fluxes    = numpy.empty((n_sets,), numpy.float64)
        syn_avg_rates     = numpy.empty((n_sets,), numpy.float64)
        
        w_n_bins = len(str(self.n_bins))
        w_n_sets = len(str(n_sets))
        
        for ibin in self.analysis_initial_bins:
            if self.bf_mode:
                trans_ibin = transdat_ds[transdat_ibin == ibin]
            else:
                trans_ibin = transdat_ds[(transdat_ibin == ibin) & transdat_in_range]
                
            for fbin in self.analysis_final_bins:
                #trans_ifbins = trans_ibin[trans_ibin['final_bin'] == fbin]
                trans_ifbins = numpy.extract(trans_ibin['final_bin'] == fbin, trans_ibin)
                dlen = len(trans_ifbins)
                
                if not dlen: continue
                
                trans_weights = trans_ifbins['final_weight']
                trans_durations = trans_ifbins['duration']
                trans_fpts      = trans_ifbins['fpt']
                trans_ibinprobs = trans_ifbins['initial_bin_pop']
                                    
                durations[ibin,fbin]['expectation'] = numpy.average(trans_durations, weights=trans_weights) * dt
                fpts[ibin,fbin]['expectation'] = numpy.average(trans_fpts, weights=trans_weights) * dt
                avg_flux = trans_weights.sum() / total_time
                fluxes[ibin,fbin]['expectation']    = avg_flux
                rates[ibin,fbin]['expectation']     = avg_flux / trans_ibinprobs.mean()
                
                for iset in range(n_sets):
                    westpa.rc.pstatus('\r  {:{w_n_bins}d}->{:<{w_n_bins}d} set {:{w_n_sets}d}/{:<{w_n_sets}d}, set size {:<20d}'\
                            .format(ibin,fbin,iset+1,n_sets,dlen, w_n_bins=w_n_bins, w_n_sets=w_n_sets), end='')
                    westpa.rc.pflush()
                    indices = numpy.random.randint(dlen, size=(dlen,))
                    #syn_weights   = trans_weights[indices]
                    #syn_durations = trans_durations[indices]
                    #syn_fpts      = trans_fpts[indices]
                    #syn_ibinprobs = trans_ibinprobs[indices]
                    syn_weights    = trans_weights.take(indices)
                    syn_durations  = trans_durations.take(indices)
                    syn_fpts       = trans_fpts.take(indices)
                    syn_ibinprobs  = trans_ibinprobs.take(indices)
                    
                    syn_avg_durations[iset] = numpy.average(syn_durations, weights=syn_weights) * dt
                    syn_avg_fpts[iset] = numpy.average(syn_fpts, weights=syn_weights) * dt
                    syn_avg_fluxes[iset] = syn_weights.sum() / total_time
                    syn_avg_rates[iset] = syn_avg_fluxes[iset] / syn_ibinprobs.mean()
                    
                    del indices, syn_weights, syn_durations, syn_ibinprobs, syn_fpts
                    
                syn_avg_durations.sort()
                syn_avg_fpts.sort()
                syn_avg_fluxes.sort()
                syn_avg_rates.sort()
                
                durations[ibin,fbin]['ci_lower'] = syn_avg_durations[lbi]
                durations[ibin,fbin]['ci_upper'] = syn_avg_durations[ubi]
                
                fpts[ibin,fbin]['ci_lower'] = syn_avg_fpts[lbi]
                fpts[ibin,fbin]['ci_upper'] = syn_avg_fpts[ubi]
                
                fluxes[ibin,fbin]['ci_lower'] = syn_avg_fluxes[lbi]
                fluxes[ibin,fbin]['ci_upper'] = syn_avg_fluxes[ubi]
                
                rates[ibin,fbin]['ci_lower'] = syn_avg_rates[lbi]
                rates[ibin,fbin]['ci_upper'] = syn_avg_rates[ubi]
                    
                del trans_weights, trans_durations, trans_ibinprobs, trans_ifbins, trans_fpts
            westpa.rc.pstatus()
            del trans_ibin
            
        for (dsname, data) in (('duration', durations), ('fpt', fpts), ('flux', fluxes), ('rate', rates)):
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
            
            if not self.bf_mode:
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
        
        if not self.bf_mode:
            self.record_data_iter_range(self.ttimes_group)
        self.record_data_binhash(self.ttimes_group)
            
    def summarize_stats(self):
        for (array, dsname, argname, title) in ((self.durations, 'duration', 'ed_stats_filename', 'event duration'),
                                                (self.fpts, 'fpt', 'fpt_stats_filename', 'first passage time'),
                                                (self.fluxes, 'flux', 'flux_stats_filename', 'flux'),
                                                (self.rates, 'rate', 'rate_stats_filename', 'rate')):
            filename = getattr(self,argname)
            if filename:
                if array is None:
                    try:
                        array = self.ttimes_group[dsname]
                    except KeyError:
                        westpa.rc.pstatus('{} data not found in {}'.format(title, self.anal_h5name))
                        continue

                self.summarize_ci(filename, array, title, self.mcbs_display_confidence,
                                  headers=(not self.suppress_headers), labels=self.print_bin_labels)
        
            
    def summarize_ci(self, filename, array, title, confidence, headers, labels):
        
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

        for ibin in self.analysis_initial_bins:
            for fbin in self.analysis_final_bins:
                mean = array[ibin,fbin]['expectation']
                lb = array[ibin,fbin]['ci_lower']
                ub = array[ibin,fbin]['ci_upper']
                ciwidth = ub - lb
                relciwidth = abs(ciwidth/mean)
                symmerr = max(mean-lb,ub-mean)
                
                outfile.write(format_2d.format(*list(map(float,(mean,lb,ub,ciwidth,relciwidth,symmerr))),
                                               ibin=ibin,fbin=fbin,mw=max_ibin_width))

    def main(self):            
        parser = argparse.ArgumentParser('w_ttimes', description=self.description)
        westpa.rc.add_args(parser)
        self.add_args(parser)
                
        args = parser.parse_args()
        westpa.rc.process_args(args, config_required = False)
        self.process_args(args)
                        
        self.check_iter_range()
        self.check_bin_selection()
        self.open_analysis_backing()
        self.ttimes_group = self.require_analysis_group('w_ttimes', replace=False)
        self.require_bin_assignments()
        self.require_transitions()
        self.gen_stats()
        self.summarize_stats()


class WTTimesWE(WTTimesBase,CommonOutputMixin,MCBSMixin,KineticsAnalysisMixin,TransitionAnalysisMixin,BinningMixin,
              IterRangeMixin,WESTDataReaderMixin,
              WESTAnalysisTool):
    description = 'Trace the WEST trajectory tree and report on transition kinetics.'
    
    def __init__(self):
        super(WTTimesWE,self).__init__()
            
                

class WTTimesBF(WTTimesBase,CommonOutputMixin,MCBSMixin,KineticsAnalysisMixin,BFTransitionAnalysisMixin,
                 BFDataManager,WESTAnalysisTool):
    description = 'Trace one or more brute force trajectories and report on transition kinetics.'
    default_chunksize = 65536*4
    
    def __init__(self):
        super(WTTimesBF,self).__init__()
        self.bf_mode = True
        self.config_required = False
        self.usecols = None
        self.input_files = None
                
    def check_iter_range(self):
        pass # do nothing, since we don't do iteration ranges for brute force
    
    def get_total_time(self):
        self.require_bf_h5file()
        return numpy.add.reduce([self.get_traj_len(traj_id)-1 for traj_id in range(self.get_n_trajs())]) * self.dt
        
