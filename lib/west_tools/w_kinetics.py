# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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

import sys, logging
from collections import deque

import numpy

import westpa
from westpa import h5io
from west.data_manager import weight_dtype
from west.data_manager import seg_id_dtype
from westpa.binning import index_dtype
from westpa.kinetics import find_macrostate_transitions
from westpa.kinetics._kinetics import _fast_transition_state_copy #@UnresolvedImport
from westpa.kinetics.matrates import estimate_rates
from westtools import (WESTSubcommand, WESTMasterCommand, WESTDataReader, IterRangeSelection,
                       ProgressIndicatorComponent)


ed_list_dtype = numpy.dtype([('istate', numpy.uint16), ('fstate', numpy.uint16), ('duration', numpy.float64),
                             ('weight', numpy.float64), ('seg_id', seg_id_dtype)])

log = logging.getLogger('westtools.w_kinetics')

class KineticsSubcommands(WESTSubcommand):
    '''Base class for common options for both kinetics schemes'''
    
    def __init__(self, parent):
        super(KineticsSubcommands,self).__init__(parent)
        self.progress = ProgressIndicatorComponent()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection() 
        self.output_file = None
        self.assignments_file = None
        
        self.do_compression = True

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
        
        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('-a', '--assignments', default='assign.h5',
                             help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                                (default: %(default)s).''')
        # default_kinetics_file will be picked up as a class attribute from the appropriate
        # subclass
        iogroup.add_argument('-o', '--output', dest='output', default=self.default_kinetics_file,
                             help='''Store results in OUTPUT (default: %(default)s).''')
        iogroup.add_argument('--no-compression', dest='compression', action='store_false',
                             help='''Do not store kinetics results compressed. This can increase disk
                             use about 100-fold, but can dramatically speed up subsequent analysis
                             for "w_kinavg matrix". Default: compress kinetics results.''')
        agroup = parser.add_argument_group('other options')
        agroup.add_argument('--config-from-file', dest='config_from_file', action='store_true', 
                            help='''Load bins/macrostates from a scheme specified in west.cfg.''')
        agroup.add_argument('--scheme-name', dest='scheme',
                            help='''Name of scheme specified in west.cfg.''')
        self.progress.add_args(parser)
        parser.set_defaults(compression=True)
        
    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)
        if args.config_from_file == False:
            self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
            self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        if args.config_from_file:
            if not args.scheme:
                raise ValueError('A scheme must be specified.')
            else:
                self.load_config_from_west(args.scheme)
        h5io.stamp_creator_data(self.output_file)
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments do not span the requested iterations')
        self.do_compression = args.compression

    def load_config_from_west(self, scheme):
        try:
            config = westpa.rc.config['west']['w_ipython']
        except:
            raise ValueError('There is no configuration file specified.')
        import os
        path = os.path.join(os.getcwd(), config['directory'], scheme)
        try:
            os.mkdir(config['directory'])
            os.mkdir(path)
        except:
            pass
        self.output_file = h5io.WESTPAH5File(os.path.join(path, 'kintrace.h5'), 'w', creating_program=True)
        self.assignments_file = h5io.WESTPAH5File(os.path.join(path, 'assign.h5'), 'r')
        

class KinTraceSubcommand(KineticsSubcommands):
    subcommand='trace'
    default_kinetics_file = 'kintrace.h5'
    help_text = 'calculate state-to-state kinetics by tracing trajectories'
    description = '''\
Calculate state-to-state rates and transition event durations by tracing
trajectories.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).

The output generated by this program is used as input for the ``w_kinavg``
tool, which converts the flux data in the output file into average rates
with confidence intervals. See ``w_kinavg trace --help`` for more 
information.

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "kintrace.h5") contains the
following datasets:

  ``/conditional_fluxes`` [iteration][state][state]
    *(Floating-point)* Macrostate-to-macrostate fluxes. These are **not**
    normalized by the population of the initial macrostate.

  ``/conditional_arrivals`` [iteration][stateA][stateB]
    *(Integer)* Number of trajectories arriving at state *stateB* in a given
    iteration, given that they departed from *stateA*.
    
  ``/total_fluxes`` [iteration][state]
    *(Floating-point)* Total flux into a given macrostate.
    
  ``/arrivals`` [iteration][state]
    *(Integer)* Number of trajectories arriving at a given state in a given
    iteration, regardless of where they originated.

  ``/duration_count`` [iteration]
    *(Integer)* The number of event durations recorded in each iteration.
    
  ``/durations`` [iteration][event duration]
    *(Structured -- see below)*  Event durations for transition events ending
    during a given iteration. These are stored as follows:
      
      istate
        *(Integer)* Initial state of transition event.
      fstate
        *(Integer)* Final state of transition event.
      duration
        *(Floating-point)* Duration of transition, in units of tau.
      weight
        *(Floating-point)* Weight of trajectory at end of transition, **not**
        normalized by initial state population.

Because state-to-state fluxes stored in this file are not normalized by
initial macrostate population, they cannot be used as rates without further
processing. The ``w_kinavg`` command is used to perform this normalization
while taking statistical fluctuation and correlation into account. See 
``w_kinavg trace --help`` for more information.  Target fluxes (total flux
into a given state) require no such normalization.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''
    

    def go(self):
        pi = self.progress.indicator
        pi.new_operation('Initializing')
        with pi:
            self.data_reader.open('r')
            nstates = self.assignments_file.attrs['nstates']
            start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
            iter_count = stop_iter - start_iter
            durations_ds = self.output_file.create_dataset('durations', 
                                                           shape=(iter_count,0), maxshape=(iter_count,None),
                                                           dtype=ed_list_dtype,
                                                           chunks=(1,15360) if self.do_compression else None,
                                                           shuffle=self.do_compression,
                                                           compression=9 if self.do_compression else None)
            durations_count_ds = self.output_file.create_dataset('duration_count',
                                                                 shape=(iter_count,), dtype=numpy.int_, shuffle=True,compression=9)
            cond_fluxes_ds = self.output_file.create_dataset('conditional_fluxes',
                                                              shape=(iter_count,nstates,nstates), dtype=weight_dtype,
                                                              chunks=(h5io.calc_chunksize((iter_count,nstates,nstates),weight_dtype)
                                                                      if self.do_compression else None),
                                                              shuffle=self.do_compression,
                                                              compression=9 if self.do_compression else None)
            total_fluxes_ds = self.output_file.create_dataset('total_fluxes',
                                                              shape=(iter_count,nstates), dtype=weight_dtype,
                                                              chunks=(h5io.calc_chunksize((iter_count,nstates),weight_dtype)
                                                                      if self.do_compression else None),
                                                              shuffle=self.do_compression,
                                                              compression=9 if self.do_compression else None)

            cond_arrival_counts_ds = self.output_file.create_dataset('conditional_arrivals',
                                                                     shape=(iter_count,nstates,nstates), dtype=numpy.uint,
                                                                     chunks=(h5io.calc_chunksize((iter_count,nstates,nstates),
                                                                                                 numpy.uint)
                                                                             if self.do_compression else None),
                                                              shuffle=self.do_compression,
                                                              compression=9 if self.do_compression else None) 
            arrival_counts_ds = self.output_file.create_dataset('arrivals',
                                                                shape=(iter_count,nstates), dtype=numpy.uint,
                                                                chunks=(h5io.calc_chunksize((iter_count,nstates),
                                                                                            numpy.uint)
                                                                        if self.do_compression else None),
                                                                shuffle=self.do_compression,
                                                                compression=9 if self.do_compression else None)

            # copy state labels for convenience
            self.output_file['state_labels'] = self.assignments_file['state_labels'][...]

            # Put nice labels on things
            for ds in (self.output_file, durations_count_ds, cond_fluxes_ds, total_fluxes_ds):
                h5io.stamp_iter_range(ds, start_iter, stop_iter)

            # Calculate instantaneous rate matrices and trace trajectories
            last_state = None
            pi.new_operation('Tracing trajectories', iter_count)
            for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
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
            

class KinMatSubcommand(KineticsSubcommands):
    subcommand='matrix'
    default_kinetics_file = 'kinmat.h5'
    help_text = 'calculate state-to-state kinetics by using the labeled bin-to-bin rate matrix'
    description = '''\
Calculate fluxes from weighted ensemble data using a matrix approach to
extrapolate to equilibrium population and flux values. This analysis is only
appropriate for simulations performed without sources and sinks.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).

The output generated by this program is used as input for the ``w_kinavg``
tool, which converts the flux data in the output file into average rates
with confidence intervals. See ``w_kinavg matrix --help`` for more 
information.

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "kinmat.h5") contains the
following dataset:

  ``/labeled_bin_fluxes`` [iteration][state][state][bin][bin]
    *(Floating-point)* Bin-to-bin (microstate) flux matrix for each trajectory
    ensemble. ``labeled_bin_fluxes[100][0][1][2][3]`` is the flux from bin
    2 to bin 3 at iteration 100, among trajectories that switched labels from
    0 to 1 within that iteration. Units are inverse tau.
                
Because fluxes stored in this file are not normalized by initial populations, 
they cannot be used as rates without further processing. The ``w_kinavg``
command is used to perform this normalization while taking statistical
fluctuation and correlation into account. See ``w_kinavg matrix --help`` for
more information.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''
    
    def __init__(self, parent):
        super(KinMatSubcommand,self).__init__(parent)
        self.window_size = None
        self.all_lags = False
    
    def add_args(self, parser):
        wgroup = parser.add_argument_group('averaging window')        
        wgroup.add_argument('-w', '--windowsize', type=int, default=1,
                            help='''Estimate kinetics over a maximum of WINDOWSIZE iterations.
                            (Default: %(default)s).''')
        wgroup.add_argument('--all-lags', action='store_true', default=False,
                            help='Use all possible lags within window of WINDOWSIZE')

    def process_args(self, args):
        self.window_size = args.windowsize
        self.all_lags = bool(args.all_lags)

    def go(self):
        pi = self.progress.indicator
        pi.new_operation('Initializing')
        with pi:
            self.data_reader.open('r')
            nbins = self.assignments_file.attrs['nbins']
            state_labels = self.assignments_file['state_labels'][...]
            state_map = self.assignments_file['state_map'][...]
            nstates = len(state_labels)
            start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
            iter_count = stop_iter - start_iter

            weights_ring = deque(maxlen=self.window_size)
            parent_ids_ring = deque(maxlen=self.window_size)
            bin_assignments_ring = deque(maxlen=self.window_size)
            label_assignments_ring = deque(maxlen=self.window_size)

            labeled_matrix_shape = (iter_count,nstates,nstates,nbins,nbins)
            unlabeled_matrix_shape = (iter_count,nbins,nbins)
            labeled_matrix_chunks = (1, nstates, nstates, nbins, nbins)
            unlabeled_matrix_chunks = (1, nbins, nbins)

            labeled_bin_fluxes_ds = self.output_file.create_dataset('labeled_bin_fluxes',
                                                                    shape=labeled_matrix_shape,
                                                                    chunks=labeled_matrix_chunks if self.do_compression else None,
                                                                    compression=9 if self.do_compression else None,
                                                                    dtype=weight_dtype)
            labeled_bin_rates_ds = self.output_file.create_dataset('labeled_bin_rates',
                                                                   shape=labeled_matrix_shape,
                                                                   chunks=labeled_matrix_chunks if self.do_compression else None,
                                                                   compression=9 if self.do_compression else None,
                                                                   dtype=weight_dtype)
            unlabeled_bin_rates_ds = self.output_file.create_dataset('bin_rates', shape=unlabeled_matrix_shape,
                                                                     chunks=unlabeled_matrix_chunks if self.do_compression else None,
                                                                     compression=9 if self.do_compression else None,
                                                                     dtype=weight_dtype)

            fluxes = numpy.empty(labeled_matrix_shape[1:], weight_dtype)
            labeled_rates = numpy.empty(labeled_matrix_shape[1:], weight_dtype)
            unlabeled_rates = numpy.empty(unlabeled_matrix_shape[1:], weight_dtype)

            for ds in (self.output_file, labeled_bin_fluxes_ds, labeled_bin_rates_ds, unlabeled_bin_rates_ds):
                h5io.stamp_iter_range(ds, start_iter, stop_iter)

            for ds in (labeled_bin_fluxes_ds, labeled_bin_rates_ds):
                h5io.label_axes(ds, ['iteration','initial state','final state','inital bin','final bin'])

            for ds in (unlabeled_bin_rates_ds,):
                h5io.label_axes(ds, ['iteration', 'initial bin', 'final bin'])

            pi.new_operation('Calculating flux matrices', iter_count)
            # Calculate instantaneous rate matrices and trace trajectories
            for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
                # Get data from the main HDF5 file
                iter_group = self.data_reader.get_iter_group(n_iter)
                seg_index = iter_group['seg_index']
                nsegs, npts = iter_group['pcoord'].shape[0:2] 
                weights = seg_index['weight']
                parent_ids = self.data_reader.parent_id_dsspec.get_iter_data(n_iter)

                # Get bin and traj. ensemble assignments from the previously-generated assignments file
                assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
                bin_assignments = numpy.require(self.assignments_file['assignments'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                                dtype=index_dtype)
                label_assignments = numpy.require(self.assignments_file['trajlabels'][assignment_iiter + numpy.s_[:nsegs,:npts]],
                                                  dtype=index_dtype)
                labeled_pops = self.assignments_file['labeled_populations'][assignment_iiter]

                # Prepare to run analysis
                weights_ring.append(weights)
                parent_ids_ring.append(parent_ids)
                bin_assignments_ring.append(bin_assignments)
                label_assignments_ring.append(label_assignments)

                # Estimate rates using bin-to-bin fluxes
                estimate_rates(nbins, state_labels,
                               weights_ring, parent_ids_ring, bin_assignments_ring, label_assignments_ring, state_map,
                               labeled_pops,
                               self.all_lags,
                               fluxes, labeled_rates, unlabeled_rates)

                # Store bin-based kinetics data
                labeled_bin_fluxes_ds[iiter] = fluxes
                labeled_bin_rates_ds[iiter] = labeled_rates
                unlabeled_bin_rates_ds[iiter] = unlabeled_rates

                # Do a little manual clean-up to prevent memory explosion
                del iter_group, weights, parent_ids, bin_assignments, label_assignments, labeled_pops
                pi.progress += 1

class WKinetics(WESTMasterCommand):
    prog = 'w_kinetics'
    subparsers_title = 'kinetics analysis schemes'
    subcommands = [KinTraceSubcommand, KinMatSubcommand]
    
    description = '''\
Perform kinetics analysis on a weighted ensemble dataset. A valid bin
assignments file including macrostate definitions must be supplied
(-a/--assignments, usually "assign.h5"); see ``w_assign --help`` for more
information.

Two schemes are available for kinetics analysis:

  trace
    Obtain kinetics information by tracing individual trajectories from 
    state to state. This scheme is always available and correct, but
    requires the global free energy landscape to be reasonably converged to
    yield correct results. (One can use, e.g., ``w_pdist evolution`` to 
    evaluate this convergence.) See ``w_kinetics trace --help`` for more
    information.
    
  matrix
    Obtain kinetics by using the bin-to-bin (microstate) rate matrix. This
    scheme is only useful for simulations without sources or sinks, and
    requires only local convergence (i.e. correct microstate rate matrix
    elements). This scheme is somewhat less sensitive to the convergence
    of the global free energy landscape, but is very sensitive to "internal
    sinks", or groups of bins that trajectories never leave once they enter.

    **This scheme is currently under development, so results may vary.**

    See ``w_kinetics matrix --help`` for more information.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

if __name__ == '__main__':
    WKinetics().main()
