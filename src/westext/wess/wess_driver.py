# Copyright (C) 2013 Joshua L. Adelman, Matthew C. Zwier, and Lillian T. Chong
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
from westext.wess.ProbAdjust import prob_adjust

EPS = numpy.finfo(numpy.float64).eps


def reduce_array(Aij):
    """Remove empty rows and columns from an array Aij and return the reduced
        array Bij and the list of non-empty states"""

    nonempty = list(range(0,Aij.shape[0]))
    eps = numpy.finfo(Aij.dtype).eps

    for i in range(0,Aij.shape[0]):
        if (Aij[i,:] < eps).all() and (Aij[:,i] < eps).all():
            nonempty.pop(nonempty.index(i))

    nne = len(nonempty)
    Bij = numpy.zeros((nne,nne))

    for i in range(0,nne):
        for j in range(0,nne):
            Bij[i,j] = Aij[nonempty[i],nonempty[j]]

    return Bij, nonempty


class WESSDriver:
    def __init__(self, sim_manager, plugin_config):
        if not sim_manager.work_manager.is_master:
            return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        self.work_manager = sim_manager.work_manager

        self.do_reweight = check_bool(plugin_config.get('do_reweighting', False))
        self.windowsize = 0.5
        self.windowtype = 'fraction'

        windowsize = plugin_config.get('window_size')
        if windowsize is not None:
            if isinstance(windowsize,float):
                self.windowsize = windowsize
                self.windowtype = 'fraction'
                if self.windowsize <= 0 or self.windowsize > 1:
                    raise ValueError('WESS parameter error -- fractional window size must be in (0,1]')
            elif isinstance(windowsize,int):
                self.windowsize = int(windowsize)
                self.windowtype = 'fixed'
            else:
                raise ValueError('WESS parameter error -- invalid window size {!r}'.format(windowsize))
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
            
        self.write_matrices = plugin_config.get('write_matrices', False)

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

        if not self.do_reweight:
            # Reweighting not requested (or not possible)
            log.debug('Reweighting not enabled')
            return

        with self.data_manager.lock:
            wess_global_group = self.data_manager.we_h5file.require_group('wess')
            last_reweighting = int(wess_global_group.attrs.get('last_reweighting',0))

        if n_iter - last_reweighting < self.reweight_period:
            # Not time to reweight yet
            log.debug('not reweighting')
            return
        else:
            log.debug('reweighting')

        mapper = we_driver.bin_mapper
        bins = we_driver.next_iter_binning
        n_bins = len(bins)
        westpa.rc.pstatus('Averaging rates')
        averager = self.get_rates(n_iter, mapper)
        binprobs = numpy.fromiter(map(operator.attrgetter('weight'),bins), dtype=numpy.float64, count=n_bins)
        orig_binprobs = binprobs.copy()

        if self.write_matrices:
            # Create storage for ourselves
            with self.data_manager.lock:
                iter_group = self.data_manager.get_iter_group(n_iter)
                try:
                    del iter_group['wess']
                except KeyError:
                    pass
    
                wess_iter_group = iter_group.create_group('wess')
                wess_iter_group.create_dataset('avg_populations', data=averager.average_populations, compression=4)
                wess_iter_group.create_dataset('unc_populations', data=averager.stderr_populations, compression=4)
                wess_iter_group.create_dataset('avg_fluxes', data=averager.average_flux, compression=4)
                wess_iter_group.create_dataset('unc_fluxes', data=averager.stderr_flux, compression=4)
                wess_iter_group.create_dataset('avg_rates', data=averager.average_rate, compression=4)
                wess_iter_group.create_dataset('unc_rates', data=averager.stderr_rate, compression=4)


        westpa.rc.pstatus('Calculating reweighting using window size of {:d}'.format(self.eff_windowsize))
        westpa.rc.pstatus('\nBin probabilities prior to reweighting:\n{!s}'.format(binprobs))
        westpa.rc.pflush()

        rij, oldindex = reduce_array(averager.average_rate)
        uij = averager.stderr_rate[numpy.ix_(oldindex,oldindex)]

        #target_regions = numpy.where(we_driver.target_state_mask)[0]
        target_regions = list(we_driver.target_states.keys())

        flat_target_regions = []
        for target_region in target_regions:
            if target_region in oldindex:  # it is possible that the target region was removed (ie if no recycling has occurred)
                new_ibin = oldindex.index(target_region)  # this is in terms of the Rij
                flat_target_regions.append(new_ibin)

        binprobs = prob_adjust(binprobs, rij, uij, oldindex, flat_target_regions)

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
                if len(bin):
                    bin.reweight(newprob)

            wess_global_group.attrs['last_reweighting'] = n_iter

        assert (abs(1 - numpy.fromiter(map(operator.attrgetter('weight'),bins), dtype=numpy.float64, count=n_bins).sum())
                        < EPS * numpy.fromiter(map(len,bins), dtype=numpy.int, count=n_bins).sum())

        westpa.rc.pflush()

