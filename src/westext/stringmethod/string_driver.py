# Copyright (C) 2013 Joshua L. Adelman
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

import numpy as np
import types

import westpa, west
from westpa import extloader
from westpa.yamlcfg import check_bool, ConfigItemMissing
from westext.stringmethod import WESTStringMethod, DefaultStringMethod
from westpa.binning import VoronoiBinMapper


class StringDriver(object):
    def __init__(self, sim_manager, plugin_config):
        super(StringDriver, self).__init__()

        if not sim_manager.work_manager.is_master:
                return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system

        # Parameters from config file
        self.windowsize = plugin_config.get('windowsize', 10)
        self.update_interval = plugin_config.get('update_interval', 10)
        self.initial_update = plugin_config.get('initial_update', 20)
        self.priority = plugin_config.get('priority', 0)

        self.write_avg_pos = check_bool(plugin_config.get('write_avgpos', True))
        self.do_update = check_bool(plugin_config.get('do_update', True))
        self.init_from_data = check_bool(plugin_config.get('init_from_data', True))

        self.dfunc = self.get_dfunc_method(plugin_config)

        # Load method to calculate average position in a bin
        # If the method is defined in an external module, correctly bind it
        ap = self.get_avgpos_method(plugin_config)
        if hasattr(ap, 'im_class'):
            self.get_avgpos = ap
        else:
            self.get_avgpos = types.MethodType(ap, self)

        # Get initial set of string centers
        centers = self.get_initial_centers()

        try:
            sm_params = self.system.sm_params
        except AttributeError as e:
            log.error('String Driver Error: system does not define sm_params. \
                        This is required and should be added to the system definition; {}'.format(e))
            raise

        # Initialize the string
        str_method = self.get_string_method(plugin_config)

        try:
            self.strings = str_method(centers, **sm_params)
        except (TypeError, AssertionError) as e:
            log.error('String Driver Error: Failed during initialization of string method: {}'.format(e))
            raise

        # Update the BinMapper
        self.update_bin_mapper()

        # Register callback
        sim_manager.register_callback(sim_manager.prepare_new_iteration, self.prepare_new_iteration, self.priority)

        westpa.rc.pstatus('-westext.stringmethod -----------------\n')
        westpa.rc.pstatus('windowsize: {}\n'.format(self.windowsize))
        westpa.rc.pstatus('update interval: {}\n'.format(self.update_interval))
        westpa.rc.pstatus('initial update: {}\n'.format(self.initial_update))
        westpa.rc.pstatus('priority: {}\n'.format(self.priority))
        westpa.rc.pstatus('write average positions: {}\n'.format(self.write_avg_pos))
        westpa.rc.pstatus('do update: {}\n'.format(self.do_update))
        westpa.rc.pstatus('initialize from WE data: {}\n'.format(self.init_from_data))
        westpa.rc.pstatus('----------------------------------------\n')
        westpa.rc.pflush()

    def dfunc(self):
        raise NotImplementedError

    def get_avgpos(self, n_iter):
        raise NotImplementedError

    def get_dfunc_method(self, plugin_config):
        try:
            methodname = plugin_config['dfunc_method']
        except KeyError:
            raise ConfigItemMissing('dfunc_method')

        dfunc_method = extloader.get_object(methodname)

        log.info('loaded stringmethod dfunc method {!r}'.format(dfunc_method))

        return dfunc_method

    def get_avgpos_method(self, plugin_config):
        try:
            methodname = plugin_config['avgpos_method']
        except KeyError:
            raise ConfigItemMissing('avgpos_method')

        if methodname.lower() == 'cartesian':
            avgpos_method = self.avgpos_cartesian
        else:
            avgpos_method = extloader.get_object(methodname)

        log.info('loaded stringmethod avgpos method {!r}'.format(avgpos_method))

        return avgpos_method

    def get_string_method(self, plugin_config):
        try:
            methodname = plugin_config['string_method']
        except KeyError:
            raise ConfigItemMissing('string_method')

        if methodname.lower() == 'default':
            str_method = DefaultStringMethod
        else:
            str_method = extloader.get_object(methodname)

        assert issubclass(str_method, WESTStringMethod)
        log.debug('loaded stringmethod string method {!r}'.format(str_method))

        return str_method

    def get_initial_centers(self):
        self.data_manager.open_backing()

        with self.data_manager.lock:
            n_iter = max(self.data_manager.current_iteration - 1, 1)
            iter_group = self.data_manager.get_iter_group(n_iter)

            # First attempt to initialize string from data rather than system
            centers = None
            if self.init_from_data:
                log.info('Attempting to initialize stringmethod from data')

                try:
                    binhash = iter_group.attrs['binhash']
                    bin_mapper = self.data_manager.get_bin_mapper(binhash)

                    centers = bin_mapper.centers

                except:
                    log.warning('Initializing string centers from data failed; Using definition in system instead.')
                    centers = self.system.bin_mapper.centers
            else:
                log.info('Initializing string centers from system definition')
                centers = self.system.bin_mapper.centers

        self.data_manager.close_backing()

        return centers

    def update_bin_mapper(self):
        '''Update the bin_mapper using the current string'''

        westpa.rc.pstatus('westext.stringmethod: Updating bin mapper\n')
        westpa.rc.pflush()

        try:
            dfargs = getattr(self.system, 'dfargs', None)
            dfkwargs = getattr(self.system, 'dfkwargs', None)
            self.system.bin_mapper = VoronoiBinMapper(self.dfunc, self.strings.centers, 
                                                      dfargs=dfargs, 
                                                      dfkwargs=dfkwargs)
        except (ValueError, TypeError) as e:
            log.error('StringDriver Error: Failed updating the bin mapper: {}'.format(e))
            raise

    def avgpos_cartesian(self, n_iter):
        '''Get average position of replicas in each bin as of n_iter for the
        the user selected update interval'''

        nbins = self.system.bin_mapper.nbins
        ndim = self.system.pcoord_ndim

        avg_pos = np.zeros((nbins, ndim), dtype=self.system.pcoord_dtype)
        sum_bin_weight = np.zeros((nbins,), dtype=self.system.pcoord_dtype)

        start_iter = max(n_iter - min(self.windowsize, n_iter), 2)
        stop_iter = n_iter + 1

        for n in range(start_iter, stop_iter):
            with self.data_manager.lock:
                iter_group = self.data_manager.get_iter_group(n)
                seg_index = iter_group['seg_index'][...]

                pcoords = iter_group['pcoord'][:,-1,:]  # Only read final point
                bin_indices = self.system.bin_mapper.assign(pcoords)
                weights = seg_index['weight']

                pcoord_w = pcoords * weights[:,np.newaxis]
                uniq_indices = np.unique(bin_indices)

                for indx in uniq_indices:
                    avg_pos[indx,:] += pcoord_w[bin_indices == indx].sum(axis=0)

                sum_bin_weight += np.bincount(bin_indices.astype(np.int), weights=weights, minlength=nbins)

        # Some bins might have zero samples so exclude to avoid divide by zero
        occ_ind = np.nonzero(sum_bin_weight)
        avg_pos[occ_ind] /= sum_bin_weight[occ_ind][:,np.newaxis]

        return avg_pos, sum_bin_weight

    def prepare_new_iteration(self):

        n_iter = self.sim_manager.n_iter

        with self.data_manager.lock:
            iter_group = self.data_manager.get_iter_group(n_iter)

            try:
                del iter_group['stringmethod']
            except KeyError:
                pass

            sm_global_group = self.data_manager.we_h5file.require_group('stringmethod')
            last_update = int(sm_global_group.attrs.get('last_update', 0))

        if n_iter - last_update < self.update_interval or n_iter < self.initial_update or not self.do_update:
            log.debug('Not updating string this iteration')
            return
        else:
            log.debug('Updating string - n_iter: {}'.format(n_iter))

        westpa.rc.pstatus('-westext.stringmethod -----------------\n')
        westpa.rc.pstatus('westext.stringmethod: Calculating average position in string images\n')
        westpa.rc.pflush()

        avg_pos, sum_bin_weight = self.get_avgpos(n_iter)

        westpa.rc.pstatus('westext.stringmethod: Updating string\n')
        westpa.rc.pflush()

        self.strings.update_string_centers(avg_pos, sum_bin_weight)

        westpa.rc.pstatus('westext.stringmethod: String lengths: {}\n'.format(self.strings.length))
        westpa.rc.pflush()

        # Update the bin definitions
        self.update_bin_mapper()

        sm_global_group.attrs['last_update'] = n_iter
