# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
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

from __future__ import print_function, division; __metaclass__ = type
import logging

# Let's suppress those numpy warnings.
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

import numpy as np
import scipy.sparse as sp

import westpa
from west.data_manager import weight_dtype, n_iter_dtype

from westpa import h5io

# From postanalysis matrix
from westpa.binning import index_dtype
from westpa.reweight import stats_process, reweight_for_c

def calc_stats(bin_assignments, weights, fluxes, populations, trans, mask, sampling_frequency, groups=None, n_groups=None):
    fluxes.fill(0.0)
    populations.fill(0.0)
    trans.fill(0)

    stats_process(bin_assignments, weights, fluxes, populations, trans, mask, interval=sampling_frequency, groups=groups, n_groups=n_groups)

class FluxMatrix():

    def w_postanalysis_matrix(self):
        pi = self.progress.indicator
        pi.new_operation('Initializing')
        
        self.data_reader.open('r')
        nbins = self.assignments_file.attrs['nbins']

        state_labels = self.assignments_file['state_labels'][...]
        state_map = self.assignments_file['state_map'][...]
        nstates = len(state_labels)

        start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter

        nfbins = nbins * nstates

        # Determine number of groups so that we can shape everything appropriately.
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            iter_group = self.data_reader.get_iter_group(n_iter)
            n_groups = 0
            for gid, seg_ids in self.generate_groups(iter_group):
                n_groups += 1

        flux_shape = (iter_count, n_groups, nfbins, nfbins)
        pop_shape = (iter_count, n_groups, nfbins)

        h5io.stamp_iter_range(self.output_file, start_iter, stop_iter)

        bin_populations_ds = self.output_file.create_dataset('bin_populations', shape=pop_shape, dtype=weight_dtype)
        h5io.stamp_iter_range(bin_populations_ds, start_iter, stop_iter)
        h5io.label_axes(bin_populations_ds, ['iteration', 'group', 'bin'])


        flux_grp = self.output_file.create_group('iterations')
        self.output_file.attrs['nrows'] = nfbins
        self.output_file.attrs['ncols'] = nfbins


        fluxes = np.empty(flux_shape[1:], weight_dtype)
        populations = np.empty(pop_shape[1:], weight_dtype)
        trans = np.empty(flux_shape[1:], np.int64)

        # Check to make sure this isn't a data set with target states
        #tstates = self.data_reader.data_manager.get_target_states(0)
        #if len(tstates) > 0:
        #    raise ValueError('Postanalysis reweighting analysis does not support WE simulation run under recycling conditions')

        pi.new_operation('Calculating flux matrices', iter_count)
        # Calculate instantaneous statistics
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            # Get data from the main HDF5 file
            iter_group = self.data_reader.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            nsegs, npts = iter_group['pcoord'].shape[0:2] 
            weights = seg_index['weight']

            # Generate groups.
            groups = np.zeros((seg_index.shape), dtype=np.long)
            n_groups = 0
            for gid, seg_ids in self.generate_groups(iter_group):
                groups[seg_ids] = gid
                n_groups += 1


            # Get bin and traj. ensemble assignments from the previously-generated assignments file
            assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
            bin_assignments = np.require(self.assignments_file['assignments'][assignment_iiter + np.s_[:nsegs,:npts]],
                                            dtype=index_dtype)

            mask_unknown = np.zeros_like(bin_assignments, dtype=np.uint16)

            macrostate_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
            macrostate_assignments = np.require(self.assignments_file['trajlabels'][macrostate_iiter + np.s_[:nsegs,:npts]],
                                        dtype=index_dtype)

            # Transform bin_assignments to take macrostate membership into account
            bin_assignments  = nstates * bin_assignments + macrostate_assignments

            mask_indx = np.where(macrostate_assignments == nstates)
            mask_unknown[mask_indx] = 1

            # Calculate bin-to-bin fluxes, bin populations and number of obs transitions
            calc_stats(bin_assignments, weights, fluxes, populations, trans, mask_unknown, self.sampling_frequency, groups, n_groups)

            # Store bin-based kinetics data
            bin_populations_ds[iiter] = populations

            flux_iter_grp = flux_grp.create_group('iter_{:08d}'.format(n_iter))
            flux_iter_grp.attrs['nrows'] = nfbins
            flux_iter_grp.attrs['ncols'] = nfbins

            # Setup sparse data structures for flux and obs
            # Moving to a 3D format has some... uh, issues, I suppose.
            # The big issue is really just going to be writing/reading.
            # Which is a large issue!  Ugh.
            fluxes_all = []
            trans_all  = []
            fluxes_m_all = []
            trans_m_all  = []
            rows_all = []
            cols_all = []
            events = 0
            events_per_group = []
            for g in range(0, fluxes.shape[0]):
                fluxes_all.append(sp.coo_matrix(fluxes[g,:]).data)
                trans_all.append(sp.coo_matrix(trans[g,:]).data)
                rows_all.append(sp.coo_matrix(fluxes[g,:]).row)
                cols_all.append(sp.coo_matrix(fluxes[g,:]).col)
                events += fluxes_all[g].shape[0]
                events_per_group.append(fluxes_all[g].shape[0])

                # assert fluxes_sp[g].nnz == trans_sp[g].nnz

            # Create a dataset filled with the appropriate number of items...
            fluxes_tw = np.zeros((events, 2))
            trans_tw = np.zeros((events,2))
            rows_tw = np.zeros((events,2))
            cols_tw = np.zeros((events,2))
            last_point = 0
            for g in range(0, fluxes.shape[0]):
                fluxes_tw[last_point:last_point+events_per_group[g],1] = g
                trans_tw[last_point:last_point+events_per_group[g],1] = g
                rows_tw[last_point:last_point+events_per_group[g],1] = g
                cols_tw[last_point:last_point+events_per_group[g],1] = g

                print(fluxes_all[g])
                fluxes_tw[last_point:last_point+events_per_group[g],0] = fluxes_all[g]
                trans_tw[last_point:last_point+events_per_group[g],0] = trans_all[g]
                rows_tw[last_point:last_point+events_per_group[g],0] = rows_all[g]
                cols_tw[last_point:last_point+events_per_group[g],0] = cols_all[g]
                last_point += events_per_group[g]

            # Now we need to handle writing the data to h5 in a hopefully fast manner.
            # It'd be nice to just... not have to worry about this stuff, but whatever.
            flux_iter_grp.create_dataset('flux'.format(g), data=fluxes_tw, dtype=weight_dtype)
            flux_iter_grp.create_dataset('obs'.format(g), data=trans_tw, dtype=np.int32)
            flux_iter_grp.create_dataset('rows'.format(g), data=rows_tw, dtype=np.int32)
            flux_iter_grp.create_dataset('cols'.format(g), data=cols_tw, dtype=np.int32)

            # Do a little manual clean-up to prevent memory explosion
            del iter_group, weights, bin_assignments, groups
            del macrostate_assignments

            pi.progress += 1

            # Check and save the number of intermediate time points; this will be used to normalize the
            # flux and kinetics to tau in w_postanalysis_reweight.
            if self.assignments_file.attrs['subsampled'] == True or self.sampling_frequency == 'iteration':
                self.output_file.attrs['npts'] = 2
            else:
                #self.output_file.attrs['npts'] = npts if self.sampling_frequency == 'timepoint' else 2
                self.output_file.attrs['npts'] = npts
