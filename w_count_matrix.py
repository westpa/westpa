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
import logging

import numpy as np
import h5py
from h5py import h5s
import sys
import traceback
import networkx as nx
import matplotlib.pyplot as plt

import westpa
from west.data_manager import weight_dtype, n_iter_dtype
from westtools import (WESTParallelTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)
from westpa import h5io
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

# These functions have been pasted in from the older tool.  They'll need to be rewritten in.
# It's just here as a template to remind me.

def __build_history_child_dictionary__(self):
    rebuild_history = 0
    for i in xrange(0, len(self.VisitedWalkers)):
        current_iter = len(self.VisitedWalkers) - i - 1
        iter = self.__format_iteration__(current_iter + 1)
        current_walker = self.VisitedWalkers[current_iter]
        walker_children = []
        seg_index = self.__return_seg_index__(iter)
        for walker in current_walker:
            child = []
            for x in xrange(0, len(seg_index)):
                try:
                    if seg_index[x]['parent_id'] == walker and (len(list(set(self.VisitedWalkers[current_iter + 1]) & set([x])))) >= 1:
                        child.append(x)
                except (KeyError):
                    pass
            try:
                self.WalkerChildren[(str(current_iter) + ":" + str(walker))] += child 
            except (KeyError):
                self.WalkerChildren[(str(current_iter) + ":" + str(walker))] = []
                self.WalkerChildren[(str(current_iter) + ":" + str(walker))] += child


class WCountMatrix(WESTParallelTool):
    prog ='w_count_matrix'
    description = '''\
            You should write something groovy here.
-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WCountMatrix, self).__init__()
        
        self.progress = ProgressIndicatorComponent()
        self.data_reader = WESTDataReader()
        self.output_filename = None
        self.iter_range = IterRangeSelection() 
        
        self.output_file = None
        self.assignments_file = None
        # Actually, I don't think we need this, truth be told.
        self.kinetics_file = None
        
        self.evolution_mode = None
        
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.progress.add_args(parser)
        self.iter_range.add_args(parser)

        iogroup = parser.add_argument_group('input/output options')

        iogroup.add_argument('-o', '--output', dest='output', default='bugger.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        iogroup.add_argument('-a', '--assignments', default='assign.h5',
                             help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                                (default: %(default)s).''')

        #cogroup = parser.add_argument_group('calculation options')
        #cogroup.add_argument('--window-frac', type=float, default=1.0,
        #                     help='''Fraction of iterations to use in each window when running in ``cumulative`` mode.
        #                     The (1 - frac) fraction of iterations will be discarded from the start of each window.''')
        
    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_filename, 'w', creating_program=True)
        self.data_reader.open('r')
        h5io.stamp_creator_data(self.output_file)

    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        with self.data_reader:
            self.iter_range.process_args(args)
        
        self.output_filename = args.output

    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            state_labels = self.assignments_file['state_labels'][...]
            state_map = self.assignments_file['state_map'][...]
            nstates = self.assignments_file.attrs['nstates']

            assert nstates == len(state_labels)

            # Taken from Josh, and from above... 
            start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
            start_pts = range(start_iter, stop_iter)


            self.successful_trajectories = {}
            pi.new_operation('Building Graph...', len(start_pts))
            # We should do all states at once, rather than force the users to run this
            # once for every state they may be interested in.  Ergo, we need to pull the
            # state assignments from the assignment file, which we have done above.
            # We're not interested in the trivial case of k = k, however.
            #for k in xrange(nstates):
            #    for j in xrange(nstates):
            #        if k != j:

                        # We're now going to begin building a state history of successful walkers, one state to state
                        # set of transitions at a time.  Probably not as fast as it can be done, but.
                        # We'll store the results in a large dictionary...
                        # Here, we want to start at the last iteration, and count until we reach iter 1.
            self.WeightGraph = nx.MultiDiGraph()
            old_parents = None
            old_children = None
            for iiter, n_iter in enumerate(xrange(stop_iter-1, start_iter-1,-1)):
                pi.progress += 1

                # Get the current iteration group and assignment iteration group...
                iter_group = self.data_reader.get_iter_group(n_iter)
                assignment_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)

                # Converting the h5io assigned iteration into the data I actually want...
                seg_index = iter_group['seg_index']
                state_assignments = self.assignments_file['trajlabels'][assignment_iiter]

                # Find all the walkers that are in the target state, and add them in to the dictionary.  We'll start adjusting this as we go and think about it.
                # We need to ignore the time information and just pull the seg ID, essentially.  I kind of don't remember which is which, however.  I think the first dimension is
                # the seg ID, and the second time.  Anyway, that's what np.where returns, here: timepoint by timepoint info.
                # But in truth, I only care that it happened sometime during this iteration, so we can ignore all that timepoint information and just return the first dimension.

                # By converting to a set and then back to a list, we remove duplicates.  There's probably a faster way with numpy, but we can worry about it later.
                # This is fast enough, for the moment, and I've marked why I'm doing it.  Anyway, we have a list of our segment IDs which are in our target states,
                # and the parents.
                # I think we're just going to use networkX to make our life a little easier here, however.
                #in_state_walkers = list(set(np.where(state_assignments == j)[0]))
                # Just return where they don't equal the non-existent state for now...
                in_state_walkers = list(set(np.where(state_assignments != nstates)[0]))
                in_state_walkers_parents = seg_index['parent_id']

                # Let's add nodes as tuples of type Iter, SegID.  We won't add any attributes, for now, although we might later.
                for i in in_state_walkers:
                    self.WeightGraph.add_node((n_iter,i))
                if old_parents != None:
                    for i in old_children:
                        self.WeightGraph.add_edge((n_iter+1,i), (n_iter, old_parents[i]))
                old_children = in_state_walkers
                old_parents = in_state_walkers_parents
            plt.figure(1,figsize=(80,80))
            nx.draw_spring(self.WeightGraph)
            plt.savefig("atlas.png",dpi=75)




if __name__ == '__main__':
    WCountMatrix().main()
