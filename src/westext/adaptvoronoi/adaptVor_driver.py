# Copyright (C) 2018 Ali Sinan Saglam
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

import westpa, west
from westpa import extloader
from westpa.yamlcfg import check_bool, ConfigItemMissing
from westpa.binning import VoronoiBinMapper


class AdaptiveVoronoiDriver:
    '''
    This plugin implements an adaptive scheme using voronoi bins from
    Zhang 2010, J Chem Phys, 132. The options exposed to the configuration
    file are:

      - av_enabled (bool, default False): Enables adaptive binning
      - max_centers (int, default 10): The maximum number of voronoi centers to be placed
      - walk_count (integer, default 5): Number of walkers per voronoi center
      - center_freq (ingeter, default 1): Frequency of center placement
      - priority (integer, default 1): Priority in the plugin order
      - dfunc_method (function, non-optional, no default): Non-optional user defined
          function that will be used to calculate distances between voronoi centers and 
          data points
      - mapper_func (function, optional): Optional user defined function for building bin
          mappers for more complicated binning schemes e.g. embedding the voronoi binning
          in a portion of the state space. If not defined the plugin will build a 
          VoronoiBinMapper with the information it has.
    '''
    def __init__(self, sim_manager, plugin_config):

        if not sim_manager.work_manager.is_master:
                return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system

        # Parameters from config file
        # this enables the adaptive voronoi, allows turning adaptive scheme off
        self.doAdaptiveVoronoi = \
             check_bool(plugin_config.get('av_enabled', False))
        # sets maximim number of centers/voronoi bins
        self.max_centers = plugin_config.get('max_centers', 10)
        # sets number of walkers per bin/voronoi center
        self.walk_count = plugin_config.get('walk_count', 5)
        # center placement frequency in number of iterations
        self.center_freq = plugin_config.get('center_freq', 1)
        # priority of the plugin (allows for order of execution)
        self.priority = plugin_config.get('priority', 0)
        # pulls the distance function that will be used by the plugin
        self.dfunc = self.get_dfunc_method(plugin_config)
        # pulls a user defined function to build the next bin mapper
        self.mapper_func = self.get_mapper_func(plugin_config)

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
            sim_manager.register_callback(sim_manager.prepare_new_iteration,
                                  self.prepare_new_iteration, self.priority)

    def dfunc(self):
        '''
        Distance function to be used by the plugin. This function
        will be used to calculate the distance between each point.
        '''
        raise NotImplementedError

    def get_dfunc_method(self, plugin_config):
        try:
            methodname = plugin_config['dfunc_method']
        except KeyError:
            raise ConfigItemMissing('dfunc_method')

        dfunc_method = extloader.get_object(methodname)

        log.info('loaded adaptive voronoi dfunc method {!r}'.format(dfunc_method))

        return dfunc_method

    def get_mapper_func(self, plugin_config):
        try:
            methodname = plugin_config['mapper_func']
        except KeyError:
            return False

        mapper_func = extloader.get_object(methodname)

        log.info('loaded adaptive voronoi mapper function {!r}'.format(mapper_func))

        return mapper_func 

    def get_initial_centers(self):
        '''
        This function pulls from the centers from either the
        previous bin mapper  or uses the definition from the
        system to calculate the number of centers
        '''
        self.data_manager.open_backing()

        with self.data_manager.lock:
            n_iter = max(self.data_manager.current_iteration - 1, 1)
            iter_group = self.data_manager.get_iter_group(n_iter)

            # First attempt to initialize voronoi centers
            # from data rather than system
            centers = None
            try:
                log.info('Voronoi centers from previous bin mapper')
                binhash = iter_group.attrs['binhash']
                bin_mapper = self.data_manager.get_bin_mapper(binhash)

                centers = bin_mapper.centers

            except:
                log.warning('Initializing voronoi centers from data failed; \
                        Using definition in system instead.')
                centers = self.system.bin_mapper.centers

        self.data_manager.close_backing()
        return centers

    def update_bin_mapper(self):
        '''Update the bin_mapper using the current set of voronoi centers'''

        westpa.rc.pstatus('westext.adaptvoronoi: Updating bin mapper\n')
        westpa.rc.pflush()

        #self.mapper_func = plugin_config.get('mapper_func', False)
        try:
            dfargs = getattr(self.system, 'dfargs', None)
            dfkwargs = getattr(self.system, 'dfkwargs', None)
            if self.mapper_func:
                # The mapper should take in 1) distance function,
                # 2) centers, 3) dfargs, 4) dfkwargs and return 
                # the mapper we want
                self.system.bin_mapper = self.mapper_func(self.dfunc, 
                                                          self.centers, 
                                                          dfargs=dfargs, 
                                                          dfkwargs=dfkwargs)
            else:
                self.system.bin_mapper = VoronoiBinMapper(self.dfunc, self.centers,
                                                          dfargs=dfargs,
                                                          dfkwargs=dfkwargs)
            self.ncenters = self.system.bin_mapper.nbins
            new_target_counts = np.empty((self.ncenters,), np.int)
            new_target_counts[...] = self.walk_count
            self.system.bin_target_counts = new_target_counts
        except (ValueError, TypeError) as e:
            log.error('AdaptiveVoronoiDriver Error: \
                    Failed updating the bin mapper: {}'.format(e))
            raise

    def update_centers(self, iter_group):
        '''
        Update the set of Voronoi centers according to
        Zhang 2010, J Chem Phys, 132. A short description
        of the algorithm can be found in the text:

        1) First reference structure is chosen randomly from
        the first set of given structure
        2) Given a set of n reference structures, for each
        configuration in the iteration the distances to each
        reference structure is calculated and the minimum
        distance is found
        3) The configuration with the minimum distance is
        selected as the next reference
        '''

        westpa.rc.pstatus('westext.adaptvoronoi: Updating Voronoi centers\n')
        westpa.rc.pflush()

        # Pull the current coordinates to find distances
        curr_pcoords = iter_group['pcoord']
        # Initialize distance array
        dists = np.zeros(curr_pcoords.shape[0])
        for iwalk, walk in enumerate(curr_pcoords):
            # Calculate distances using the provided function
            # and find the distance to the closest center
            dists[iwalk] = min(self.dfunc(walk[-1], self.centers))
        # Find the maximum of the minimum distances
        max_ind = np.where(dists == dists.max())
        # Use the maximum progress coordinate as our next center
        self.centers = np.vstack((self.centers, 
             curr_pcoords[max_ind[0][0]][-1]))

    def prepare_new_iteration(self):

        n_iter = self.sim_manager.n_iter

        with self.data_manager.lock:
            iter_group = self.data_manager.get_iter_group(n_iter)

        # Check if we are at the correct frequency for updating the bin mapper
        if n_iter % self.center_freq == 0:
            # Check if we still need to add more centers
            if self.ncenters < self.max_centers:
                # First find the center to add
                self.update_centers(iter_group)
                # Update the bin mapper with the new center
                self.update_bin_mapper()
