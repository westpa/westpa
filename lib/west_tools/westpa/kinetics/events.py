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

def _group_eval_block(**future_kwargs):
    # Our rate estimator is a little more complex, so we've defined a custom evaluation block for it,
    # instead of just using the block evalutors that we've imported.
    results = []
    '''for istate in xrange(nstates):
        for jstate in xrange(nstates):
            if istate == jstate: continue
            kwargs = { 'istate' : istate, 'jstate': jstate }
            # Why are we sending in the total population dataset, instead of a sliced one?
            # It's a requirement of our estimator; we need to pull from any given i to j state in order to properly normalize
            # and avoid i to j rate constants which are affected by a third state k.
            # That is, we need the populations for both i and j, and it's easier to just send in the entire dataset.
            dataset = {'dataset': data_input['dataset'][:, istate, jstate], 'pops': data_input['pops'] }
            ci_res = mcbs_ci_correl(dataset,estimator=sequence_macro_flux_to_rate,
                                    alpha=mcbs_alpha,n_sets=mcbs_nsets,autocorrel_alpha=mcbs_acalpha,
                                    subsample=numpy.mean, do_correl=do_correl, mcbs_enable=mcbs_enable, estimator_kwargs=kwargs)
            results.append((name, iblock, istate, jstate, (start,stop) + ci_res))'''
    nstates = future_kwargs['nstates']
    iiter = future_kwargs['iiter']
    parent_ids = future_kwargs['parent_ids']
    last_state = future_kwargs['last_state']
    state = _fast_transition_state_copy(iiter, nstates, parent_ids, last_state)
    del(future_kwargs['iiter'])
    del(future_kwargs['parent_ids'])
    del(future_kwargs['last_state'])

    # This is almost certainly not fast.
    cond_fluxes = numpy.zeros((nstates,nstates), weight_dtype)
    total_fluxes = numpy.zeros((nstates,), weight_dtype)
    cond_counts = numpy.zeros((nstates,nstates), numpy.uint)
    total_counts = numpy.zeros((nstates,), numpy.uint)
    durations = []
    future_kwargs.update(dict(macro_fluxes=cond_fluxes,
                         macro_counts=cond_counts,
                         target_fluxes=total_fluxes,
                         target_counts=total_counts,
                         state=state,
                         durations=durations))
    gid = future_kwargs['gid']
    del(future_kwargs['gid'])
    find_macrostate_transitions(**future_kwargs)
    results.append((gid, cond_fluxes, total_fluxes, cond_counts, total_counts, durations, state))
    return results

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
        # Determine the number of groups!
        n_groups = 0
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            # Get data from the main HDF5 file to determine the number of groups.
            # Hopefully, this is a fast operation.
            iter_group = self.data_reader.get_iter_group(n_iter)
            cn_groups = self.number_of_groups(iter_group)
            n_groups = max(n_groups, cn_groups[0])
            if cn_groups == n_groups:
                group_ids = cn_groups[1]

        durations_ds = self.output_file.replace_dataset('durations', 
                                                       shape=(iter_count,n_groups,0), maxshape=(iter_count,n_groups,None),
                                                       dtype=ed_list_dtype,
                                                       chunks=(1,n_groups,15360) if self.do_compression else None,
                                                       shuffle=self.do_compression,
                                                       compression=9 if self.do_compression else None)
        durations_count_ds = self.output_file.replace_dataset('duration_count',
                                                             shape=(iter_count,n_groups,), dtype=numpy.int_, shuffle=True,compression=9)
        cond_fluxes_ds = self.output_file.replace_dataset('conditional_fluxes',
                                                          shape=(iter_count,n_groups,nstates,nstates), dtype=weight_dtype,
                                                          chunks=(h5io.calc_chunksize((iter_count,n_groups,nstates,nstates),weight_dtype)
                                                                  if self.do_compression else None),
                                                          shuffle=self.do_compression,
                                                          compression=9 if self.do_compression else None)
        total_fluxes_ds = self.output_file.replace_dataset('total_fluxes',
                                                          shape=(iter_count,n_groups,nstates), dtype=weight_dtype,
                                                          chunks=(h5io.calc_chunksize((iter_count,n_groups,nstates),weight_dtype)
                                                                  if self.do_compression else None),
                                                          shuffle=self.do_compression,
                                                          compression=9 if self.do_compression else None)

        cond_arrival_counts_ds = self.output_file.replace_dataset('conditional_arrivals',
                                                                 shape=(iter_count,n_groups,nstates,nstates), dtype=numpy.uint,
                                                                 chunks=(h5io.calc_chunksize((iter_count,n_groups,nstates,nstates),
                                                                                             numpy.uint)
                                                                         if self.do_compression else None),
                                                          shuffle=self.do_compression,
                                                          compression=9 if self.do_compression else None) 
        arrival_counts_ds = self.output_file.replace_dataset('arrivals',
                                                            shape=(iter_count,n_groups,nstates), dtype=numpy.uint,
                                                            chunks=(h5io.calc_chunksize((iter_count,n_groups,nstates),
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
        pi.new_operation('Tracing trajectories', iter_count*n_groups)
        last_state = {}
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            iter_group = self.data_reader.get_iter_group(n_iter)
            nsegs, npts = iter_group['pcoord'].shape[0:2] 
            futures = []
            assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
            bin_assignments = numpy.require(self.assignments_file['assignments'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                            dtype=index_dtype)
            label_assignments = numpy.require(self.assignments_file['trajlabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                              dtype=index_dtype)
            state_assignments = numpy.require(self.assignments_file['statelabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                              dtype=index_dtype)
            
            for gid, seg_ids in self.generate_groups(iter_group):
                if gid not in last_state:
                    last_state[gid] = None
                # Get data from the main HDF5 file
                seg_ids = list(seg_ids)
                seg_index = iter_group['seg_index'][seg_ids]
                weights = seg_index['weight']
                #parent_ids = seg_index['parent_id']
                parent_ids = self.data_reader.parent_id_dsspec.get_iter_data(n_iter)[seg_ids]
                
                # Get bin and traj. ensemble assignments from the previously-generated assignments file
                # Prepare to run analysis
                '''cond_fluxes = numpy.zeros((nstates,nstates), weight_dtype)
                total_fluxes = numpy.zeros((nstates,), weight_dtype)
                cond_counts = numpy.zeros((nstates,nstates), numpy.uint)
                total_counts = numpy.zeros((nstates,), numpy.uint)
                durations = []'''

                # Estimate macrostate fluxes and calculate event durations using trajectory tracing
                # state is opaque to the find_macrostate_transitions function            
                # We'll need to sort through this a bit differently, as we're changing it from per iteration to per iteration AND group.
                # As it is, if I just leave it here in the loop, it'll error out.
                # We probably just need to copy and store the last state per GROUP.

                future_kwargs = dict(nstates=nstates, weights=weights,
                                     label_assignments=label_assignments[seg_ids],
                                     state_assignments=state_assignments[seg_ids],
                                     dt=1.0/(npts-1),
                                     iiter=iiter,
                                     parent_ids=parent_ids,
                                     last_state=last_state[gid],
                                     gid=gid)
                
                futures.append(self.work_manager.submit(_group_eval_block, kwargs=future_kwargs))

                #find_macrostate_transitions(nstates, weights, label_assignments, state_assignments, 1.0/(npts-1), state,
                #                            cond_fluxes, cond_counts, total_fluxes, total_counts, durations)
                
            for future in self.work_manager.as_completed(futures):
                results = future.get_result(discard=True)[0]
                gid, cond_fluxes, total_fluxes, cond_counts, total_counts, durations, state = results
                last_state[gid] = state
                # Store trace-based kinetics data
                cond_fluxes_ds[iiter,gid] = cond_fluxes
                total_fluxes_ds[iiter,gid] = total_fluxes
                arrival_counts_ds[iiter,gid] = total_counts
                cond_arrival_counts_ds[iiter,gid] = cond_counts
                
                durations_count_ds[iiter,gid] = len(durations)
                if len(durations) > 0:
                    durations_ds.resize((iter_count,n_groups, max(len(durations), durations_ds.shape[1])))
                    durations_ds[iiter,gid,:len(durations)] = durations
                        
                # Do a little manual clean-up to prevent memory explosion
                #del weights, parent_ids, bin_assignments, label_assignments, state, cond_fluxes, total_fluxes
                pi.progress += 1
            del iter_group
