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
                weights = iter_group['seg_index']['weight']
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
                    # Walker, then timepoint?  For the state assignment
                    self.WeightGraph.add_node((n_iter,i), weight=weights[i], state_assignments = list(set(state_assignments[i,:])), seg_id=i, iteration=n_iter)
                if old_parents != None:
                    for i in old_children:
                        #This is correct.  The iteration is meant to indicate forward progress.
                        self.WeightGraph.add_edge((n_iter+1,i), (n_iter, old_parents[i]), iteration=n_iter)
                old_children = in_state_walkers
                old_parents = in_state_walkers_parents
            # Some silly code to test to see whether I built the code correctly or not.  Neat stuff, though...
            #plt.figure(1,figsize=(80,80))
            #nx.draw_spring(self.WeightGraph)
            #plt.savefig("test.pdf",dpi=750)
            
            # We're now going to make copies of the original graph, and prune them out for trajectories that don't begin and end in state k and j, respectively.
            # This should technically re-implement the functionality of w_kinetics, if it's done correctly, but that's what we need.
            pi.new_operation('Building the state by state graphs...', len(start_pts))
            self.StateGraphs = {}
            for k in xrange(nstates):
                for j in xrange(nstates):
                    if k != j:
                        pruned_nodes = 0
                        k_nodes = []
                        self.StateGraphs[k,j] = self.WeightGraph.copy()
                        for iiter, n_iter in enumerate(xrange(stop_iter-1, start_iter-1,-1)):
                            pi.progress += 1
                            iter_nodes = self.StateGraphs[k,j].nodes(data=True)
                            niter_nodes = []
                            for i in iter_nodes:
                                if i[1]['iteration'] == n_iter:
                                    niter_nodes.append(i)
                            #if n_iter == stop_iter-1:
                            for i in niter_nodes:
                                if k in i[1]['state_assignments']:
                                    k_nodes.append(i)
                            # We've added the nodes that are in the appropriate state to an ever growing node list, and will continue to do so.
                            # Now, we check other nodes against the state list; are they in the appropriate state list?
                            # And if not, is there a path from them to one of the appropriate k_nodes?
                            # If not, remove them.
                            for i in niter_nodes:
                                if k not in i[1]['state_assignments']: 
                                    for z in k_nodes:
                                        if z[0] != i[0]:
                                            try:
                                                print(z[0], i[0])
                                                path = nx.bidirectional_dijkstra(self.StateGraphs[k,j],i[0],z[0])
                                                print(path)
                                            except(nx.exception.NetworkXNoPath):
                                                # If we hit here, we're done and need to stop comparing.
                                                print("To remove a node is a great sin.")
                                                pruned_nodes += 1
                                                self.StateGraphs[k,j].remove_node(i[0])
                                                break

                        # We shouldn't be removing more nodes than we ever had to begin with,
                        assert (pruned_nodes + len(self.StateGraphs[k,j].nodes())) == len(self.WeightGraph.nodes())
'''
            pi.new_operation('Iterating over the graph...', len(start_pts))
            n_paths = 0
            for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
                pi.progress += 1
                # We'll only return the edges from the current iteration (or only work with them, at any rate)
                # We can do this by specifying nbunch to the edges and asserting that iteration == niter,
                # or by sorting through the quite frankly not so nice output.  Hm.
                # Specifying nbunch is not likely to be viable, given that we'll be trimming these graphs for nodes not in successful pathways
                # before we get to this state.
                # We'll determine a better way to do this later.
                iter_edges = self.WeightGraph.edges(data=True)
                niter_edges = []
                for i in iter_edges:
                    if i[2] == {'iteration': n_iter}:
                        niter_edges.append(i)
                print("This is the new iteration! " + str(n_iter))
                if n_iter == 1:
                    n_paths = len(niter_edges)
                else:
                    if niter_edges != [] and old_niter_edges != []:
                        n_paths += len(niter_edges) - len(old_niter_edges)
                old_niter_edges = niter_edges
                print(len(niter_edges), len(old_niter_edges))
                print(niter_edges, old_niter_edges)
                print(n_paths)
'''



if __name__ == '__main__':
    WCountMatrix().main()
