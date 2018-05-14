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

from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)

import numpy as np

import westpa, west
from westpa import extloader
from westpa.yamlcfg import check_bool, ConfigItemMissing
from westpa.binning import VoronoiBinMapper


class AdaptiveVoronoiDriver:
    def __init__(self, sim_manager, plugin_config):

        if not sim_manager.work_manager.is_master:
                return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        
        # Parameters from config file
        self.doAdaptiveVoronoi = check_bool(plugin_config.get('av_enabled', False))
        self.max_centers = plugin_config.get('max_centers', 10)
        self.walk_count = plugin_config.get('walk_count', 5)
        self.center_freq = plugin_config.get('center_freq', 1)
        self.priority = plugin_config.get('priority', 0)

        self.dfunc = self.get_dfunc_method(plugin_config)

        # Get initial set of Voronoi centers
        self.centers = self.get_initial_centers()
        self.ncenters = len(self.centers)

        # Update the BinMapper
        self.update_bin_mapper()

        westpa.rc.pstatus('-adaptive voronoi mapping --------------\n')
        westpa.rc.pstatus('enabled: {}\n'.format(self.doAdaptiveVoronoi))
        westpa.rc.pstatus('max centers: {}\n'.format(self.max_centers))
        westpa.rc.pstatus('center adding freq: {}\n'.format(self.center_freq))
        westpa.rc.pstatus('centers: {}\n'.format(self.centers))
        westpa.rc.pstatus('----------------------------------------\n')
        westpa.rc.pflush()

        # Register callback
        if self.doAdaptiveVoronoi:
            sim_manager.register_callback(sim_manager.prepare_new_iteration, self.prepare_new_iteration, self.priority)

    def dfunc(self):
        raise NotImplementedError

    def get_dfunc_method(self, plugin_config):
        try:
            methodname = plugin_config['dfunc_method']
        except KeyError:
            raise ConfigItemMissing('dfunc_method')

        dfunc_method = extloader.get_object(methodname)

        log.info('loaded adaptive voronoi dfunc method {!r}'.format(dfunc_method))

        return dfunc_method

    def get_initial_centers(self):
        self.data_manager.open_backing()

        with self.data_manager.lock:
            n_iter = max(self.data_manager.current_iteration - 1, 1)
            iter_group = self.data_manager.get_iter_group(n_iter)

            # First attempt to initialize string from data rather than system
            centers = None
            try:
                log.info('Voronoi centers from previous bin mapper')
                binhash = iter_group.attrs['binhash']
                bin_mapper = self.data_manager.get_bin_mapper(binhash)

                centers = bin_mapper.centers

            except:
                log.warning('Initializing string centers from data failed; Using definition in system instead.')
                centers = self.system.bin_mapper.centers

        self.data_manager.close_backing()
        return centers

    def update_bin_mapper(self):
        '''Update the bin_mapper using the current set of voronoi centers'''

        westpa.rc.pstatus('westext.adaptvoronoi: Updating bin mapper\n')
        westpa.rc.pflush()


        try:
            dfargs = getattr(self.system, 'dfargs', None)
            dfkwargs = getattr(self.system, 'dfkwargs', None)
            self.system.bin_mapper = VoronoiBinMapper(self.dfunc, self.centers, 
                                                      dfargs=dfargs, 
                                                      dfkwargs=dfkwargs)
            new_target_counts = np.empty((self.system.bin_mapper.nbins,), np.int)
            new_target_counts[...] = self.walk_count
            self.system.bin_target_counts = new_target_counts
            self.ncenters = self.system.bin_mapper.nbins
        except (ValueError, TypeError) as e:
            log.error('AdaptiveVoronoiDriver Error: Failed updating the bin mapper: {}'.format(e))
            raise

    def update_centers(self, iter_group):
        '''Update the set of Voronoi centers'''

        westpa.rc.pstatus('westext.adaptvoronoi: Updating Voronoi centers\n')
        westpa.rc.pflush()

        curr_pcoords = iter_group['pcoord']
        dists = np.zeros(curr_pcoords.shape[0])
        for iwalk, walk in enumerate(curr_pcoords):
            dists[iwalk] = min(self.dfunc(walk[-1], self.centers))
        max_ind = np.where(dists==dists.max())
        #print("prev centers")
        #print(self.centers, self.centers.shape)
        #print(curr_pcoords[max_ind[0][0]][-1])
        self.centers = np.vstack((self.centers, curr_pcoords[max_ind[0][0]][-1]))
        #print("new centers")
        #print(self.centers, self.centers.shape)

    def prepare_new_iteration(self):

        n_iter = self.sim_manager.n_iter

        with self.data_manager.lock:
            iter_group = self.data_manager.get_iter_group(n_iter)

        if n_iter % self.center_freq == 0:
            if self.ncenters < self.max_centers:
                self.update_centers(iter_group)
                # Update the bin definitions
                self.update_bin_mapper()
