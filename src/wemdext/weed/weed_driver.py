from __future__ import division; __metaclass__ = type

import logging
log = logging.getLogger(__name__)

import math,numpy
from itertools import izip

from wemd.util.miscfn import vgetattr
import wemd
from wemdext.weed.ProbAdjustEquil import probAdjustEquil

class WEEDDriver:
    def __init__(self, sim_manager):
        if sim_manager.work_manager.mode != 'master':
            return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        
        self.do_reweight = wemd.rc.config.get_bool('weed.do_equilibrium_reweighting', False)
        self.windowsize = wemd.rc.config.get_int('weed.window_size', 0)
        self.reweight_period = wemd.rc.config.get_int('weed.reweight_period', 0)
        
        if sim_manager.system.target_states and self.do_reweight:
            log.warning('equilibrium reweighting requested but target states (sinks) present; reweighting disabled')
            self.do_reweight = False 
        else:
            sim_manager.register_callback(sim_manager.prepare_new_segments, self.prepare_new_segments)    
        
            
    def get_rates(self, n_iter, bins):
        '''Get rates and associated uncertainties as of n_iter, according to the window size the user
        has selected (self.windowsize)'''
        
        n_bins = len(bins)
        
        eff_windowsize = min(n_iter,self.windowsize or n_iter)
        rates = numpy.empty((eff_windowsize,n_bins,n_bins), numpy.float64)
                
        n_used = 0
        for n in xrange(n_iter, n_iter-eff_windowsize, -1):
            log.debug('considering iteration {:d}'.format(n))
            iter_group = self.data_manager.get_iter_group(n)
            rates_ds = iter_group['bin_rates']
            if rates_ds.shape != rates.shape[1:]:
                # A bin topology change means we can't go any farther back
                self.sim_manager.status_stream.write(('Rate matrix for iteration {:d} is of the wrong shape; '
                                                      +'stopping accumulation of average rate data.\n').format(n))
                break
            
            rates[n_used] = rates_ds[...]
            n_used += 1

        avg_rates = rates[:n_used].mean(axis=0)
        if n_used == 1:
            unc_rates = avg_rates.copy()
        else:
            unc_rates = rates[:n_used].std(axis=0) / math.sqrt(n_used)
        return (avg_rates, unc_rates, n_used)
        
        
    def prepare_new_segments(self, n_iter, old_segments, new_segments):
        status_stream = self.sim_manager.status_stream
        
        if not self.do_reweight:
            # Reweighting not requested (or not possible)
            log.debug('equilibrium reweighting not enabled') 
            return

        iter_group = self.data_manager.get_iter_group(n_iter)
        
        region_set = self.sim_manager.system.region_set
        region_set.clear()
        bins = region_set.get_all_bins()
        n_bins = len(bins)
        initial_pcoords = [segment.pcoord[0] for segment in new_segments]
        for (bin, segment) in izip(region_set.map_to_bins(initial_pcoords), new_segments):
            bin.add(segment)
        
        try:
            del iter_group['weed']
        except KeyError:
            pass
        
        weed_iter_group = iter_group.create_group('weed')
        avg_rates_ds = weed_iter_group.create_dataset('avg_rates', shape=(n_bins,n_bins), dtype=numpy.float64)
        unc_rates_ds = weed_iter_group.create_dataset('unc_rates', shape=(n_bins,n_bins), dtype=numpy.float64)
        weed_global_group = self.data_manager.h5file.require_group('weed')
        last_reweighting = long(weed_global_group.attrs.get('last_reweighting', 0))
        
        if n_iter - last_reweighting < self.reweight_period:
            # Not time to reweight yet
            log.debug('not reweighting')
            return
        else:
            log.debug('reweighting')
        
        avg_rates, unc_rates, eff_windowsize = self.get_rates(n_iter, bins)
        avg_rates_ds[...] = avg_rates
        unc_rates_ds[...] = unc_rates
        
        binprobs = iter_group['bin_populations'][-1,:]
        assert numpy.allclose(binprobs, vgetattr('weight', bins, numpy.float64))
        orig_binprobs = binprobs.copy()
        
        status_stream.write('Calculating equilibrium reweighting using window size of {:d}\n'.format(eff_windowsize))
        status_stream.write('\nBin probabilities prior to reweighting:\n{!s}\n'.format(binprobs))
        self.sim_manager.flush_status()
        
        probAdjustEquil(binprobs, avg_rates, unc_rates)
        
        # Check to see if reweighting has set non-zero bins to zero probability (should never happen)
        assert (~((orig_binprobs > 0) & (binprobs == 0))).all(), 'populated bin reweighted to zero probability'
        
        # Check to see if reweighting has set zero bins to nonzero probability (may happen)
        z2nz_mask = (orig_binprobs == 0) & (binprobs > 0) 
        if (z2nz_mask).any():
            status_stream.write('Reweighting would assign nonzero probability to an empty bin; not reweighting this iteration.\n')
            status_stream.write('Empty bins assigned nonzero probability: {!s}.\n'
                                .format(numpy.array_str(numpy.arange(n_bins)[z2nz_mask])))
        else:
            status_stream.write('\nBin populations after reweighting:\n{!s}\n'.format(binprobs))
            for (bin, newprob) in izip(bins, binprobs):
                bin.reweight(newprob)
            weed_global_group.attrs['last_reweighting'] = n_iter
                
        region_set.clear()
            
        assert abs(1 - vgetattr('weight', new_segments, numpy.float64).sum()) < 1.0e-15*len(new_segments)
