'''
w_fluxanl: Obtain rate constants from flux into target states.

Monte Carlo bootstrapping is used to obtain confidence intervals on the
resulting rates.
'''

from __future__ import division, print_function

import argparse
import numpy
import scipy.signal
from itertools import izip
import westpa, oldtools

from oldtools.aframe import WESTAnalysisTool,WESTDataReaderMixin,IterRangeMixin,MCBSMixin

import logging
log = logging.getLogger('w_fluxanl')

rstat_dtype = numpy.dtype([('count', numpy.uint),
                           ('flux', numpy.float64),
                           ])

class WFluxanl(MCBSMixin,IterRangeMixin,WESTDataReaderMixin,WESTAnalysisTool):
    def __init__(self):
        super(WFluxanl,self).__init__()
        
        self.tau = None
        self.suppress_headers = None
        self.wfl_group = None
        
    def collect_flux_data(self):
        '''Assemble flux data from each iteration into a few data sets'''
        
        n_iters = self.last_iter - self.first_iter + 1
        n_targets = self.get_iter_group(self.first_iter)['recycling'].shape[0]
        
        westpa.rc.pstatus('Collecting fluxes and counts for {:d} target state(s).'.format(n_targets))
        rstats = numpy.empty((n_iters,n_targets), dtype=rstat_dtype)

        for itarget in xrange(n_targets):
            for (iiter,n_iter) in enumerate(xrange(self.first_iter,self.last_iter+1)):
                recycling = self.get_iter_group(n_iter)['recycling'][itarget]
                rstats[iiter,itarget] = (recycling['count'], recycling['flux'])
                
        rstats['flux'] /= self.tau
        rstats_ds = self.wfl_group.create_dataset('arrivals', data=rstats, compression='gzip')
        rstats_ds.attrs['tau'] = self.tau
        rstats_ds.attrs['first_iter'] = self.first_iter
        rstats_ds.attrs['iter_step'] = 1
        
        self.wfl_group.attrs['n_targets'] = n_targets
        
    def calc_autocorrel(self):
        '''Calculate the autocorrelation function of flux values and report the estimated
        time required to reach steady state.'''
        
        westpa.rc.pstatus('Calculating flux autocorrelation time...')
        n_targets = self.wfl_group.attrs['n_targets']
        lbi, ubi = self.calc_ci_bound_indices()
        dlen = self.wfl_group['arrivals'].shape[0]
        
        acorr_ds = self.wfl_group.create_dataset('autocorrel', shape=(dlen,n_targets,3), dtype=numpy.float64)
        
        for itarget in xrange(n_targets):
            fluxes = self.wfl_group['count'][:,itarget]['flux']
            ffluxes = fluxes - fluxes.mean()
            
            acorr = scipy.signal.correlate(ffluxes, ffluxes)[-dlen:]
            acorr /= acorr.max()
            acorr_bounds = numpy.empty((dlen,2), numpy.float64)
            syn_acorr = numpy.empty((self.mcbs_nsets, len(acorr)), numpy.float64)
            
            for iset in xrange(self.mcbs_nsets):
                westpa.rc.pstatus('\r  Set {:d}/{:d}'.format(iset+1,self.mcbs_nsets), end='')
                indices = numpy.random.randint(dlen, size=(dlen,))
                syn_fluxes = ffluxes[indices]
                syn_acorr[iset,:] = scipy.signal.correlate(syn_fluxes,syn_fluxes)[-dlen:]
                syn_acorr[iset,:] /= syn_acorr[iset,:].max()
                westpa.rc.pflush()
                
            westpa.rc.pstatus()
            for ilag in xrange(1,dlen):
                syn_acorr[:,ilag].sort()
                acorr_bounds[ilag,0] = syn_acorr[lbi,ilag]
                acorr_bounds[ilag,1] = syn_acorr[ubi,ilag]
                                
            acorr_ds[:,itarget,0] = acorr
            acorr_ds[:,itarget,1:] = acorr_bounds
            
            
        
    def calc_blocked_flux_cis(self):
        '''Calculate confidence intervals of average flux within blocks of iterations'''
        
        n_blocks = self.n_iter_blocks()
        n_targets = self.wfl_group.attrs['n_targets']
        lbi, ubi = self.calc_ci_bound_indices()
        
        block_bounds = numpy.empty((n_blocks,2), numpy.min_scalar_type(self.last_iter))
        flux_cis = numpy.empty((n_blocks,n_targets), dtype=oldtools.aframe.mcbs.ciinfo_dtype)
        
        all_fluxes = self.wfl_group['arrivals']['flux']
        
        westpa.rc.pstatus('Calculating flux confidence intervals...')
        for iblock, (blk_begin, blk_end) in enumerate(self.iter_block_iter()):
            westpa.rc.pstatus('\r  Iterations [{:d},{:d})'.format(blk_begin, blk_end), end='')
            #print("Averaging over iterations [{:d},{:d}).".format(blk_begin,blk_end))
            block_bounds[iblock] = blk_begin,blk_end-1
            iibegin = blk_begin - self.first_iter
            iiend   = blk_end - self.first_iter
            
            
            fluxes = all_fluxes[iibegin:iiend]
            dlen = len(fluxes)
            
            syn_avg_flux = numpy.empty((self.mcbs_nsets,n_targets), numpy.float64)

            for iset in xrange(self.mcbs_nsets):
                indices = numpy.random.randint(dlen, size=(dlen,))
                syn_fluxes = fluxes[indices]
                syn_avg_flux[iset] = syn_fluxes.mean(axis=0)
                
            flux_cis[iblock]['expectation'] = fluxes.mean(axis=0)
            for itarget in xrange(n_targets):
                syn_avg_flux[:,itarget].sort()
                flux_cis[iblock,itarget]['ci_lower'] = syn_avg_flux[lbi,itarget]
                flux_cis[iblock,itarget]['ci_upper'] = syn_avg_flux[ubi,itarget]

            westpa.rc.pflush()
        
        self.wfl_group['blocked_fluxes'] = flux_cis
        self.wfl_group['blocked_iter_bounds'] = block_bounds
        
        self.wfl_group['blocked_fluxes'].attrs['mcbs_alpha'] = self.mcbs_alpha
        self.wfl_group['blocked_fluxes'].attrs['mcsb_nsets'] = self.mcbs_nsets
        
        westpa.rc.pstatus()
        
    def calc_cumul_flux_cis(self):
        '''Calculate confidence intervals of average flux as a simulation progresses'''
        
        n_blocks = self.n_iter_blocks()
        n_targets = self.wfl_group.attrs['n_targets']
        lbi, ubi = self.calc_ci_bound_indices()
        
        block_bounds = numpy.empty((n_blocks,2), numpy.min_scalar_type(self.last_iter))
        block_bounds[:,0] = self.first_iter
        flux_cis = numpy.empty((n_blocks,n_targets), dtype=oldtools.aframe.mcbs.ciinfo_dtype)
        
        all_fluxes = self.wfl_group['arrivals']['flux']
        
        westpa.rc.pstatus('Calculating cumulative flux confidence intervals...')
        for iblock, (blk_begin, blk_end) in enumerate(self.iter_block_iter()):
            westpa.rc.pstatus('\r  Iterations [{:d},{:d})'.format(blk_begin, blk_end), end='')
            block_bounds[iblock] = self.first_iter,blk_end-1
            iiend   = blk_end - self.first_iter
            fluxes = all_fluxes[:iiend]
            dlen = len(fluxes)
            
            syn_avg_flux = numpy.empty((self.mcbs_nsets,n_targets), numpy.float64)

            for iset in xrange(self.mcbs_nsets):
                indices = numpy.random.randint(dlen, size=(dlen,))
                syn_fluxes = fluxes[indices]
                syn_avg_flux[iset] = syn_fluxes.mean(axis=0)
                
            flux_cis[iblock]['expectation'] = fluxes.mean(axis=0)
            for itarget in xrange(n_targets):
                syn_avg_flux[:,itarget].sort()
                flux_cis[iblock,itarget]['ci_lower'] = syn_avg_flux[lbi,itarget]
                flux_cis[iblock,itarget]['ci_upper'] = syn_avg_flux[ubi,itarget]

            westpa.rc.pflush()
        
        self.wfl_group['cumul_fluxes'] = flux_cis
        self.wfl_group['cumul_iter_bounds'] = block_bounds
        
        self.wfl_group['cumul_fluxes'].attrs['mcbs_alpha'] = self.mcbs_alpha
        self.wfl_group['cumul_fluxes'].attrs['mcsb_nsets'] = self.mcbs_nsets
        
        westpa.rc.pstatus()
        
        
            
    def write_output(self, dataset, output_pattern):
        if not output_pattern:
            return

        n_targets = self.wfl_group.attrs['n_targets']
        iter_bounds = self.wfl_group['{}_iter_bounds'.format(dataset)][...]
        mw = len(str(iter_bounds.max()))
        fmt = '    '.join(['{first_iter:{mw}d}', '{last_iter:{mw}d}'] + ['{:<24.16g}']*6) + '\n'
        for itarget in xrange(n_targets):
            stats_output = file(output_pattern % itarget, 'wt')
            if not args.suppress_headers:
                stats_output.write('''\
# Weighted ensemble average flux
# Target index:      {itarget:d}
# tau:               {tau:g}
# Confidence level:  {confidence:s}%
# ----
# column 0: first iteration of averaging window
# column 1: last iteration of averaging window
# column 2: average flux
# column 3: lower bound of CI
# column 4: upper bound of CI
# column 5: width of CI
# column 6: relative width of CI [abs(width/average)]
# column 7: symmetrized error [max(upper bound - average, average - lower bound)]
# ----            
'''.format(itarget=int(itarget), tau=float(self.tau), confidence=self.mcbs_display_confidence))
            
            flux_cis = self.wfl_group['{}_fluxes'.format(dataset)][:,itarget]
            assert len(flux_cis) == len(iter_bounds)
            
            for ((first_iter, last_iter), flux) in izip(iter_bounds, flux_cis):
                avg = flux['expectation']
                ci_lower = flux['ci_lower']
                ci_upper = flux['ci_upper']
                ci_width = ci_upper - ci_lower
                rel_ci_width = abs(ci_width/avg) if avg else 0
                symmerr = max(ci_upper-avg,avg-ci_lower)                
                stats_output.write(fmt.format(*map(float,map(float,[avg,ci_lower,ci_upper,ci_width,rel_ci_width,symmerr])),
                                              first_iter=long(first_iter),last_iter=long(last_iter),mw=mw))            
            stats_output.close()
                
wfl = WFluxanl()

parser = argparse.ArgumentParser('w_fluxanl')
westpa.rc.add_args(parser)
wfl.add_args(parser)

cgroup = parser.add_argument_group('calculation options')
cgroup.add_argument('-t', '--tau', dest='tau', type=float, default=1.0,
                    help='Scale rates to correspond to a propagation/resampling timestep (default: 1.0).')
cgroup.add_argument('--calc-autocorrel', dest='calc_autocorrel', action='store_true',
                    help='Calculate flux autocorrelation function.')

ogroup = parser.add_argument_group('output options')
ogroup.add_argument('--output-blocked', dest='output_pattern_blocked', default='avgflux_%s_blocked.txt', metavar='PATTERN',
                    help='Store output for block averaged fluxes in PATTERN, which must contain a printf-style pattern '
                        +'which will be replaced by the name/index of the target state to which results correspond '
                        +' (default: %(default)s).')
ogroup.add_argument('--output-cumul', dest='output_pattern_cumul', default='avgflux_%s_cumul.txt', metavar='PATTERN',
                    help='Store output for cumulative average flux in PATTERN, which must contain a printf-style pattern '
                        +'which will be replaced by the name/index of the target state to which results correspond '
                        +' (default: %(default)s).')
ogroup.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')

args = parser.parse_args()
westpa.rc.process_args(args, config_required=False)
wfl.process_args(args)
wfl.tau = args.tau
wfl.suppress_headers = args.suppress_headers

wfl.check_iter_range()
wfl.open_analysis_backing()
wfl.wfl_group = wfl.require_analysis_group('w_fluxanl', replace=True)

wfl.collect_flux_data()
if args.calc_autocorrel:
    wfl.calc_autocorrel()
wfl.calc_blocked_flux_cis()
wfl.calc_cumul_flux_cis()
wfl.write_output('blocked', args.output_pattern_blocked)
wfl.write_output('cumul', args.output_pattern_cumul)

