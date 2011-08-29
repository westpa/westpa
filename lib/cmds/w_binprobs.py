from __future__ import division, print_function
import argparse, numpy, warnings

import logging
log = logging.getLogger('w_binprobs')

import wemd

from wemdtools.aframe import WEMDAnalysisTool, BinningMixin, DataReaderMixin, IterRangeMixin, MCBSMixin
from wemdtools.aframe.mcbs import calc_ci_bound_indices

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

# Upcalls move left to right
# If some complex dependencies exist, one can always override the process_args function and
# call parent classes' process_args() manually in an order that makes sense.
class WBinprobs(MCBSMixin, BinningMixin, IterRangeMixin, DataReaderMixin, WEMDAnalysisTool):
    def calc_average_binprobs(self):
        '''Calculate average bin populations over blocks of iterations, with MCBS error bars.'''
        
        lbi, ubi = calc_ci_bound_indices(self.mcbs_nsets, self.mcbs_alpha)
        if self.iter_step == 1:
            wemd.rc.pstatus('Calculating per-iteration average bin populations...')
        else:
            wemd.rc.pstatus('Calculating average bin populations in blocks of {:d} iterations...'.format(self.iter_step))
        all_pops = self.binning_h5group['bin_populations'][...]
        pcoord_len = self.get_pcoord_len(self.first_iter)
        syn_avg_pops = numpy.empty((self.mcbs_nsets, self.n_bins), numpy.float64)

        first_n_iters = xrange(self.first_iter, 
                               self.last_iter+self.iter_step if self.last_iter % self.iter_step else self.last_iter, 
                               self.iter_step)

        populations = numpy.empty((len(first_n_iters), self.n_bins), dtype=ciinfo_dtype)
        iter_bounds = numpy.empty((len(first_n_iters), 2), numpy.min_scalar_type(self.last_iter))
        
        # Iterate over blocks
        mw = len(str(self.last_iter))
        for iblock, n_iter in enumerate(first_n_iters):
            wemd.rc.pstatus('\r  Iterations {:{mw}d} -- {:<{mw}d}'.format(n_iter,n_iter+self.iter_step-1,mw=mw), end='')
            
            iifirst = n_iter-self.first_iter
            iilast  = iifirst+self.iter_step
            
            iter_pops = all_pops[iifirst:iilast]
            n_iters_block = len(iter_pops)        # number of iterations in this block
            
            iter_bounds[iblock] = (n_iter, n_iter+n_iters_block-1)

            poplen = (len(iter_pops)-1)*(pcoord_len-1) + pcoord_len
            blockpops = numpy.empty((poplen,self.n_bins), numpy.float64)
            
            # Pool per-timepoint, per-bin populations without double counting the time points at
            # iteration boundaries
            for iiter in xrange(n_iters_block-1):
                blockpops[iiter*(pcoord_len-1):(iiter+1)*(pcoord_len-1),:] = iter_pops[iiter,:-1,:]
            blockpops[(n_iters_block-1)*(pcoord_len-1):,:] = iter_pops[n_iters_block-1]
            dlen = len(blockpops)
            
            populations[iblock]['expectation'] = blockpops.mean(axis=0) # per-bin average populations for this block

            for iset in xrange(self.mcbs_nsets):
                indices = numpy.random.randint(dlen,size=(dlen,))
                syn_pops = blockpops[indices]
                syn_avg_pops[iset, :] = syn_pops.mean(axis=0)
            
            for ibin in xrange(self.n_bins):
                syn_avg_pops[:,ibin].sort()
            
            populations[iblock]['ci_lower'] = syn_avg_pops[lbi,:]
            populations[iblock]['ci_upper'] = syn_avg_pops[ubi,:]
            
            wemd.rc.pflush()
            
        group = self.anal_h5file['w_binprobs']
        group.attrs['mcbs_alpha'] = self.mcbs_alpha
        group.attrs['mcbs_nsets'] = self.mcbs_nsets
        
        # Write a numerically-friendly matrix of bin probabilities to the HDF5 file
        group['average_binprobs'] = populations
        group.create_dataset('iter_bounds', data=iter_bounds)
        dsattrs = group['average_binprobs'].attrs
        dsattrs['first_iter'] = self.first_iter
        dsattrs['iter_step'] = self.iter_step
        
        wemd.rc.pstatus()

        # Explicit deletes so that memory is returned to the system immediately
        # instead of at the next garbage collection            
        del populations, iter_bounds, all_pops, iter_pops, blockpops, syn_avg_pops
        
    def write_output(self):
        if not self.output_filename:
            return
        
        output_file = open(self.output_filename, 'wt')
        
        if not self.suppress_headers:
            output_file.write('''\
# Average WEMD bin probabilities
# Iterations {first_iter} -- {last_iter} (inclusive)
# Confidence level: {confidence}%
# Number of bootstrap data sets: {mcbs_nsets}
# ----
'''.format(first_iter=self.first_iter, last_iter=self.last_iter, confidence=self.mcbs_display_confidence,
           mcbs_nsets = self.mcbs_nsets))
            
            if self.do_write_bin_labels:
                self.write_bin_labels(output_file)
                output_file.write('----\n')
                
            output_file.write('''\
# column 0: first iteration of averaging window
# column 1: last iteration of averaging window
# column 2: bin index
# column 3: average population
# column 4: lower bound of confidence interval
# column 5: upper bound of confidence interval
# column 6: width of confidence interval
# column 7: relative width of confidence interval [abs(width/average)]
# column 8: symmetrized error bar [max(upper bound - average, average - lower bound)]
# ----
''')
        
        
        group = self.anal_h5file['w_binprobs']
        iter_bounds = group['iter_bounds'][...]
        populations = group['average_binprobs'][...]
        last_iter = iter_bounds[-1,-1]
        
        iw = len(str(last_iter))
        bw = len(str(self.n_bins-1))
        fmt = '    '.join(['{first_iter:{iw}d}',
                           '{last_iter:{iw}d}',
                           '{ibin:{bw}d}',
                           '{avg:< 22.15g}',
                           '{ci_lower:< 22.15g}',
                           '{ci_upper:< 22.15g}',
                           '{ci_width:< 22.15g}',
                           '{rel_ci_width:< 22.15g}',
                           '{sym_err:< 22.15g}\n'])
        
        for iblock in xrange(len(iter_bounds)):
            first_iter, last_iter = iter_bounds[iblock]
            for ibin in xrange(self.n_bins):
                avg, ci_lower, ci_upper = populations[iblock, ibin]
                ci_width = ci_upper - ci_lower
                rel_ci_width = abs(ci_width / avg) if avg != 0 else 0
                sym_err = max(ci_upper-avg, avg-ci_lower)
                output_file.write(fmt.format(iw=iw, bw=bw,
                                             first_iter=long(first_iter), last_iter=long(last_iter),
                                             ibin=long(ibin),
                                             avg=float(avg), ci_lower=float(ci_lower), ci_upper=float(ci_upper),
                                             ci_width=float(ci_width), rel_ci_width=float(rel_ci_width), sym_err=float(sym_err)))
        
        output_file.close()
        
wbp = WBinprobs()

parser = argparse.ArgumentParser('w_binprobs')
wemd.rc.add_args(parser)
wbp.add_args(parser)

ogroup = parser.add_argument_group('output options')
ogroup.add_argument('-o', '--output', default='binprobs.txt', 
                    help='''Write average bin probabilities to OUTPUT (default: %(default)s).''')
ogroup.add_argument('--no-headers', dest='suppress_headers', action='store_true',
                    help='''Do not write (commented) headers to output files.''')
ogroup.add_argument('--bin-labels', dest='write_bin_labels', action='store_true',
                    help='''Write bin labels to output files.''')

args = parser.parse_args()

wemd.rc.process_args(args, config_required=False)
wbp.process_args(args)

wbp.output_filename = args.output
wbp.suppress_headers = bool(args.suppress_headers)
wbp.do_write_bin_labels = bool(args.write_bin_labels)

wbp.check_iter_range()
wbp.open_analysis_backing()
wbp.require_bin_assignments()
wbp.require_analysis_group('w_binprobs', replace=True)
wbp.calc_average_binprobs()
wbp.write_output()


    