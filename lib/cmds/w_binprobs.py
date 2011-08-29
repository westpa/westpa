from __future__ import division, print_function
import os, sys, argparse, numpy
from math import ceil, log10

import logging
log = logging.getLogger('w_binprobs')

import wemd, wemdtools

from wemdtools.aframe import WEMDAnalysisTool, BinningMixin, DataReaderMixin, IterRangeMixin, MCBSMixin
from wemdtools.aframe.mcbs import calc_ci_bound_indices

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

# Upcalls move left to right
# If some complex dependencies exist, one can always override the process_common_args function and
# call parent classes' process_common_args() manually in an order that makes sense.
class WBinprobs(MCBSMixin, BinningMixin, IterRangeMixin, DataReaderMixin, WEMDAnalysisTool):
        def calc_average_binprobs(self):
        '''Calculate average bin populations over blocks of iterations, with MCBS error bars.'''
        
        lbi, ubi = calc_ci_bound_indices(self.mcbs_nsets, self.mcbs_alpha)
        if self.iter_step == 1:
            wemd.rc.pstatus('Calculating per-iteration average bin populations')
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
        wemd.rc.pstatus('Calculating average bin populations...')
        for iblock, n_iter in enumerate(first_n_iters):
            wemd.rc.pstatus('\r  Iteration {:d}'.format(n_iter), end='')
            
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
            
        group = self.require_analysis_group('w_binprobs')
        group['average_binprobs'] = populations
        group.create_dataset('average_binprobs_iters', data=iter_bounds, compression='gzip') 
        
        attrs = group['average_binprobs'].attrs
        attrs['mcbs_alpha'] = self.mcbs_alpha
        attrs['mcbs_nsets'] = self.mcbs_nsets 
        wemd.rc.pstatus()

        # Explicit deletes so that memory is returned to the system immediately
        # instead of at the next garbage collection            
        del populations, iter_bounds, all_pops, iter_pops, blockpops, syn_avg_pops
        
wbp = WBinprobs()

parser = argparse.ArgumentParser('w_binprobs')
wemd.rc.add_common_args(parser)
wbp.add_common_args(parser)

args = parser.parse_args()

wemd.rc.process_common_args(args, config_required=False)
wbp.process_common_args(args)
wbp.check_iter_range()
wbp.open_analysis_backing()
wbp.require_bin_assignments()

wbp.require_analysis_group('w_binprobs', replace=True)
wbp.calc_average_binprobs()

    