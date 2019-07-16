# Copyright (C) 2013 Matthew C. Zwier, Joshua L. Adelman and Lillian T. Chong
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


import logging
log = logging.getLogger(__name__)

import numpy
import operator


import westpa, west
from westpa.yamlcfg import check_bool
from westpa.kinetics import RateAverager
from westext.weed.ProbAdjustEquil import probAdjustEquil

EPS = numpy.finfo(numpy.float64).eps


class WEEDDriver:
    def __init__(self, sim_manager, plugin_config):
        if not sim_manager.work_manager.is_master:
            return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        self.work_manager = sim_manager.work_manager

        self.do_reweight = check_bool(plugin_config.get('do_equilibrium_reweighting', False)
                           or plugin_config.get('do_reweighting', False))
        self.windowsize = 0.5
        self.windowtype = 'fraction'

        windowsize = plugin_config.get('window_size')
        if windowsize is not None:
            if isinstance(windowsize, float):
                self.windowsize = windowsize
                self.windowtype = 'fraction'
                if self.windowsize <= 0 or self.windowsize > 1:
                    raise ValueError('WEED parameter error -- fractional window size must be in (0,1]')
            elif isinstance(windowsize,int):
                self.windowsize = int(windowsize)
                self.windowtype = 'fixed'
            else:
                raise ValueError('WEED parameter error -- invalid window size {!r}'.format(windowsize))
        log.info('using window size of {!r} ({})'.format(self.windowsize, self.windowtype))

        self.max_windowsize = plugin_config.get('max_window_size')
        if self.max_windowsize is not None:
            log.info('Using max windowsize of {:d}'.format(self.max_windowsize))

        self.reweight_period = plugin_config.get('reweight_period', 0)
        self.priority = plugin_config.get('priority', 0)

        self.rate_calc_queue_size = plugin_config.get('rate_calc_queue_size', 1)
        self.rate_calc_n_blocks = plugin_config.get('rate_calc_n_blocks', 1)

        if self.do_reweight:
            sim_manager.register_callback(sim_manager.prepare_new_iteration,self.prepare_new_iteration, self.priority)

    def get_rates(self, n_iter, mapper):
        '''Get rates and associated uncertainties as of n_iter, according to the window size the user
        has selected (self.windowsize)'''

        if self.windowtype == 'fraction':
            if self.max_windowsize is not None:
                eff_windowsize = min(self.max_windowsize,int(n_iter * self.windowsize))
            else:
                eff_windowsize = int(n_iter * self.windowsize)
        else: # self.windowtype == 'fixed':
            eff_windowsize = min(n_iter, self.windowsize or 0)

        averager = RateAverager(mapper, self.system, self.data_manager, self.work_manager)
        averager.calculate(max(1, n_iter-eff_windowsize), n_iter+1, self.rate_calc_n_blocks, self.rate_calc_queue_size)
        self.eff_windowsize = eff_windowsize

        return averager

    def prepare_new_iteration(self):
        n_iter = self.sim_manager.n_iter
        we_driver = self.sim_manager.we_driver
        
        if we_driver.target_states and self.do_reweight:
            log.warning('equilibrium reweighting requested but target states (sinks) present; reweighting disabled')
            return 

        if not self.do_reweight:
            # Reweighting not requested
            log.debug('equilibrium reweighting not enabled') 
            return

        with self.data_manager.lock:
            weed_global_group = self.data_manager.we_h5file.require_group('weed')
            last_reweighting = int(weed_global_group.attrs.get('last_reweighting', 0))

        if n_iter - last_reweighting < self.reweight_period:
            # Not time to reweight yet
            log.debug('not reweighting')
            return
        else:
            log.debug('reweighting')

        mapper = we_driver.bin_mapper
        bins = we_driver.next_iter_binning
        n_bins = len(bins)

        # Create storage for ourselves
        with self.data_manager.lock:
            iter_group = self.data_manager.get_iter_group(n_iter)
            try:
                del iter_group['weed']
            except KeyError:
                pass

            weed_iter_group = iter_group.create_group('weed')
            avg_populations_ds = weed_iter_group.create_dataset('avg_populations', shape=(n_bins,), dtype=numpy.float64)
            unc_populations_ds = weed_iter_group.create_dataset('unc_populations', shape=(n_bins,), dtype=numpy.float64)
            avg_flux_ds = weed_iter_group.create_dataset('avg_fluxes', shape=(n_bins,n_bins), dtype=numpy.float64)
            unc_flux_ds = weed_iter_group.create_dataset('unc_fluxes', shape=(n_bins,n_bins), dtype=numpy.float64)
            avg_rates_ds = weed_iter_group.create_dataset('avg_rates', shape=(n_bins,n_bins), dtype=numpy.float64)
            unc_rates_ds = weed_iter_group.create_dataset('unc_rates', shape=(n_bins,n_bins), dtype=numpy.float64)

        averager = self.get_rates(n_iter, mapper)

        with self.data_manager.flushing_lock():
            avg_populations_ds[...] = averager.average_populations
            unc_populations_ds[...] = averager.stderr_populations
            avg_flux_ds[...] = averager.average_flux
            unc_flux_ds[...] = averager.stderr_flux
            avg_rates_ds[...] = averager.average_rate
            unc_rates_ds[...] = averager.stderr_rate

            binprobs = numpy.fromiter(map(operator.attrgetter('weight'),bins), dtype=numpy.float64, count=n_bins)
            orig_binprobs = binprobs.copy()

        westpa.rc.pstatus('Calculating equilibrium reweighting using window size of {:d}'.format(self.eff_windowsize))
        westpa.rc.pstatus('\nBin probabilities prior to reweighting:\n{!s}'.format(binprobs))
        westpa.rc.pflush()

        probAdjustEquil(binprobs, averager.average_rate, averager.stderr_rate)
        
        # Check to see if reweighting has set non-zero bins to zero probability (should never happen)
        assert (~((orig_binprobs > 0) & (binprobs == 0))).all(), 'populated bin reweighted to zero probability'
        
        # Check to see if reweighting has set zero bins to nonzero probability (may happen)
        z2nz_mask = (orig_binprobs == 0) & (binprobs > 0) 
        if (z2nz_mask).any():
            westpa.rc.pstatus('Reweighting would assign nonzero probability to an empty bin; not reweighting this iteration.')
            westpa.rc.pstatus('Empty bins assigned nonzero probability: {!s}.'
                                .format(numpy.array_str(numpy.arange(n_bins)[z2nz_mask])))
        else:
            westpa.rc.pstatus('\nBin populations after reweighting:\n{!s}'.format(binprobs))
            for (bin, newprob) in zip(bins, binprobs):
                bin.reweight(newprob)

            weed_global_group.attrs['last_reweighting'] = n_iter

        assert (abs(1 - numpy.fromiter(map(operator.attrgetter('weight'),bins), dtype=numpy.float64, count=n_bins).sum())
                < EPS * numpy.fromiter(map(len,bins), dtype=numpy.int, count=n_bins).sum())

        westpa.rc.pflush()
