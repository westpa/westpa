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


# Let's suppress those numpy warnings.
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

import numpy

from west.data_manager import weight_dtype, n_iter_dtype, seg_id_dtype
from westpa import h5io

# From w_kinetics.
from westtools.dtypes import ed_list_dtype
from westpa.binning import index_dtype
from westpa.kinetics._kinetics import _fast_transition_state_copy #@UnresolvedImport
from westpa.kinetics import find_macrostate_transitions

# The old w_kinetics
class WKinetics():
    def w_kinetics(self):
        pi = self.progress.indicator
        pi.new_operation('Initializing')

        self.data_reader.open('r')
        self.open_files()
        nstates = self.assignments_file.attrs['nstates']
        start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter
        durations_ds = self.output_file.replace_dataset('durations', 
                                                       shape=(iter_count,0), maxshape=(iter_count,None),
                                                       dtype=ed_list_dtype,
                                                       chunks=(1,15360) if self.do_compression else None,
                                                       shuffle=self.do_compression,
                                                       compression=9 if self.do_compression else None)
        durations_count_ds = self.output_file.replace_dataset('duration_count',
                                                             shape=(iter_count,), dtype=numpy.int_, shuffle=True,compression=9)
        cond_fluxes_ds = self.output_file.replace_dataset('conditional_fluxes',
                                                          shape=(iter_count,nstates,nstates), dtype=weight_dtype,
                                                          chunks=(h5io.calc_chunksize((iter_count,nstates,nstates),weight_dtype)
                                                                  if self.do_compression else None),
                                                          shuffle=self.do_compression,
                                                          compression=9 if self.do_compression else None)
        total_fluxes_ds = self.output_file.replace_dataset('total_fluxes',
                                                          shape=(iter_count,nstates), dtype=weight_dtype,
                                                          chunks=(h5io.calc_chunksize((iter_count,nstates),weight_dtype)
                                                                  if self.do_compression else None),
                                                          shuffle=self.do_compression,
                                                          compression=9 if self.do_compression else None)

        cond_arrival_counts_ds = self.output_file.replace_dataset('conditional_arrivals',
                                                                 shape=(iter_count,nstates,nstates), dtype=numpy.uint,
                                                                 chunks=(h5io.calc_chunksize((iter_count,nstates,nstates),
                                                                                             numpy.uint)
                                                                         if self.do_compression else None),
                                                          shuffle=self.do_compression,
                                                          compression=9 if self.do_compression else None) 
        arrival_counts_ds = self.output_file.replace_dataset('arrivals',
                                                            shape=(iter_count,nstates), dtype=numpy.uint,
                                                            chunks=(h5io.calc_chunksize((iter_count,nstates),
                                                                                        numpy.uint)
                                                                    if self.do_compression else None),
                                                            shuffle=self.do_compression,
                                                            compression=9 if self.do_compression else None)

        # copy state labels for convenience
        self.output_file.replace_dataset('state_labels', data=self.assignments_file['state_labels'][...])

        # Put nice labels on things
        for ds in (self.output_file, durations_count_ds, cond_fluxes_ds, total_fluxes_ds):
            h5io.stamp_iter_range(ds, start_iter, stop_iter)

        # Calculate instantaneous rate matrices and trace trajectories
        last_state = None
        pi.new_operation('Tracing trajectories', iter_count)
        for iiter, n_iter in enumerate(range(start_iter, stop_iter)):
            # Get data from the main HDF5 file
            iter_group = self.data_reader.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            nsegs, npts = iter_group['pcoord'].shape[0:2] 
            weights = seg_index['weight']
            #parent_ids = seg_index['parent_id']
            parent_ids = self.data_reader.parent_id_dsspec.get_iter_data(n_iter)
            
            # Get bin and traj. ensemble assignments from the previously-generated assignments file
            assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
            bin_assignments = numpy.require(self.assignments_file['assignments'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                            dtype=index_dtype)
            label_assignments = numpy.require(self.assignments_file['trajlabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                              dtype=index_dtype)
            state_assignments = numpy.require(self.assignments_file['statelabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                              dtype=index_dtype)
            
            # Prepare to run analysis
            cond_fluxes = numpy.zeros((nstates,nstates), weight_dtype)
            total_fluxes = numpy.zeros((nstates,), weight_dtype)
            cond_counts = numpy.zeros((nstates,nstates), numpy.uint)
            total_counts = numpy.zeros((nstates,), numpy.uint)
            durations = []

            # Estimate macrostate fluxes and calculate event durations using trajectory tracing
            # state is opaque to the find_macrostate_transitions function            
            state = _fast_transition_state_copy(iiter, nstates, parent_ids, last_state)
            find_macrostate_transitions(nstates, weights, label_assignments, state_assignments, 1.0/(npts-1), state,
                                        cond_fluxes, cond_counts, total_fluxes, total_counts, durations)
            last_state = state
            
            # Store trace-based kinetics data
            cond_fluxes_ds[iiter] = cond_fluxes
            total_fluxes_ds[iiter] = total_fluxes
            arrival_counts_ds[iiter] = total_counts
            cond_arrival_counts_ds[iiter] = cond_counts
            
            durations_count_ds[iiter] = len(durations)
            if len(durations) > 0:
                durations_ds.resize((iter_count, max(len(durations), durations_ds.shape[1])))
                durations_ds[iiter,:len(durations)] = durations
                    
            # Do a little manual clean-up to prevent memory explosion
            del iter_group, weights, parent_ids, bin_assignments, label_assignments, state, cond_fluxes, total_fluxes
            pi.progress += 1
