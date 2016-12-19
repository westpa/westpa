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
import scipy as sp
from scipy import signal
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

# For now, we'll borrow from Matt and see how this works out...
from mclib import mcbs_ci,autocorrel_elem
from matplotlib.pyplot import plot
from matplotlib.mlab import cohere


def _assign_label_pop(n_iter, lb, ub, mapper, nstates, state_map, last_labels, parent_id_dsspec, weight_dsspec, pcoord_dsspec):    

    nbins = len(state_map)-1
    parent_ids = parent_id_dsspec.get_iter_data(n_iter,index_exp[lb:ub])
    weights = weight_dsspec.get_iter_data(n_iter,index_exp[lb:ub])
    pcoords = pcoord_dsspec.get_iter_data(n_iter,index_exp[lb:ub])
    
    assignments, trajlabels = assign_and_label(lb, ub, parent_ids,
                                               mapper.assign, nstates, state_map, last_labels, pcoords)
    pops = numpy.zeros((nstates+1,nbins+1), weight_dtype)
    accumulate_labeled_populations(weights, assignments, trajlabels, pops)
    return (assignments, trajlabels, pops, lb, ub)


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
 
    def build_graph_iteration(self, n_iter, nstates, nbins, state_map, last_labels):
        ''' Method to encapsulate the segment slicing (into n_worker slices) and parallel job submission
            Submits job(s), waits on completion, splices them back together
            Returns: assignments, trajlabels, pops for this iteration'''

        futures = []

        iter_group = self.data_reader.get_iter_group(n_iter)
        nsegs, npts = iter_group['pcoord'].shape[:2]
        n_workers = self.work_manager.n_workers or 1
        assignments = numpy.empty((nsegs, npts), dtype=index_dtype)
        trajlabels = numpy.empty((nsegs, npts), dtype=index_dtype)
        pops = numpy.zeros((nstates+1,nbins+1), dtype=weight_dtype)

        #Submit jobs to work manager
        blocksize = nsegs // n_workers
        if nsegs % n_workers > 0:
            blocksize += 1

        def task_gen():
            if __debug__:
                checkset = set()
            for lb in xrange(0, nsegs, blocksize):
                ub = min(nsegs, lb+blocksize)
                if __debug__:
                    checkset.update(set(xrange(lb,ub)))
                args = ()
                kwargs = dict(n_iter=n_iter,
                              lb=lb, ub=ub, mapper=self.binning.mapper, nstates=nstates, state_map=state_map,
                              last_labels=last_labels, 
                              parent_id_dsspec=self.data_reader.parent_id_dsspec, 
                              weight_dsspec=self.data_reader.weight_dsspec,
                              pcoord_dsspec=self.dssynth.dsspec)
                yield (_assign_label_pop, args, kwargs)

                #futures.append(self.work_manager.submit(_assign_label_pop, 
                #kwargs=)
            if __debug__:
                assert checkset == set(xrange(nsegs)), 'segments missing: {}'.format(set(xrange(nsegs)) - checkset)

        #for future in self.work_manager.as_completed(futures):
        for future in self.work_manager.submit_as_completed(task_gen(), queue_size=self.max_queue_len):
            assign_slice, traj_slice, slice_pops, lb, ub = future.get_result(discard=True)
            assignments[lb:ub, :] = assign_slice
            trajlabels[lb:ub, :] = traj_slice
            pops += slice_pops
            del assign_slice, traj_slice, slice_pops

        del futures
        return (assignments, trajlabels, pops)

    def prune_graph(self, graph, k, k_nodes, n_iter, pruned_nodes):
        iter_nodes = graph.nodes(data=True)
        niter_nodes = []
        for i in iter_nodes:
            #print(i)
            try:
                if i[1]['iteration'] == n_iter:
                    niter_nodes.append(i)
            except:
                #print(i)
                pass
        for i in niter_nodes:
            #print(i[1]['state_assignments'])
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
                            path = nx.bidirectional_dijkstra(graph,i[0],z[0])
                            # We need to break if this is successful, as well, because a single path to ANY point is good enough to leave the node in,
                            # Otherwise, we tend to do things like a cut node out if it happens to not go to a particular successful end point.
                            break
                        except(nx.exception.NetworkXNoPath):
                            # If we hit here, we're done and need to stop comparing.
                            #print("Not in network.")
                            pruned_nodes += 1
                            graph.remove_node(i[0])
                            break
                if k_nodes == []:
                    #print("Cutting some bitches!")
                    graph.remove_node(i[0])
                    pruned_nodes += 1
        return graph, k_nodes, pruned_nodes

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
                npcoord = iter_group['pcoord']
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
                # Just return where they don't equal the non-existent state for now... also, the unknown state.
                in_state_walkers = list(set(np.where(state_assignments < nstates)[0]))
                in_state_walkers_parents = seg_index['parent_id']

                # Let's add nodes as tuples of type Iter, SegID.  We won't add any attributes, for now, although we might later.
                for i in in_state_walkers:
                    # Walker, then timepoint?  For the state assignment
                    self.WeightGraph.add_node((n_iter,i), weight=weights[i], state_assignments=list(set(state_assignments[i,:])), seg_id=i, iteration=n_iter, pcoord=npcoord[i,:])
                if old_parents != None:
                    for i in old_children:
                        # This is correct.  The iteration is meant to indicate forward progress.
                        #self.WeightGraph.add_edge((n_iter+1,i), (n_iter, old_parents[i]), iteration=n_iter)
                        # Actually, I think I had it backwards.  n_iter was alright, though.  I might have to go back and check all the other routines for correctness now, however.  But the direction is definitely good!
                        # Oh, these are edges, not nodes.
                        self.WeightGraph.add_edge((n_iter, old_parents[i]), (n_iter+1,i), iteration=n_iter)
                old_children = in_state_walkers
                old_parents = in_state_walkers_parents
            # Some silly code to test to see whether I built the code correctly or not.  Neat stuff, though...
            #plt.figure(1,figsize=(80,80))
            #nx.draw_spring(self.WeightGraph)
            #plt.savefig("test.pdf",dpi=750)
            
            # We're now going to make copies of the original graph, and prune them out for trajectories that don't begin and end in state k and j, respectively.
            # This should technically re-implement some of the idea of w_kinetics (i.e., pathfinding), if it's done correctly, but that's what we need.

            pi.new_operation('Building the state by state graphs...', len(start_pts))
            self.StateGraphs = {}
            pruned_nodes = {}
            for k in xrange(nstates):
                for j in xrange(nstates):
                    if k != j:
                        pruned_nodes[k,j] = 0
                        k_nodes = []
                        self.StateGraphs[k,j] = self.WeightGraph.copy()
                        for iiter, n_iter in enumerate(xrange(stop_iter-1, start_iter-1,-1)):
                            print(n_iter)
                            pi.progress += 1

                            # I can do this better, I think.
                            print("Run 1 " + str(k) + " " + str(j))
                            self.StateGraphs[k,j], k_nodes, pruned_nodes[k,j] = self.prune_graph(self.StateGraphs[k,j], j, k_nodes, n_iter, pruned_nodes[k,j])
                            #print(k_nodes)

                        assert (pruned_nodes[k,j] + len(self.StateGraphs[k,j].nodes())) == len(self.WeightGraph.nodes())

            # ... anyway, do it again!  Just turn it around, and do it for states ending in the j state, this time.  I may have flipped this around, but I think that's alright.
            for k in xrange(nstates):
                for j in xrange(nstates):
                    if k != j:
                        k_nodes = []
                        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
                            print(n_iter)
                            pi.progress += 1

                            # I can do this better, I think.
                            print("Run 2 " + str(k) + " " + str(j))
                            self.StateGraphs[k,j], k_nodes, pruned_nodes[k,j] = self.prune_graph(self.StateGraphs[k,j], k, k_nodes, n_iter, pruned_nodes[k,j])
                            #print(k_nodes)

                        print(len(self.StateGraphs[k,j]))
                        assert (pruned_nodes[k,j] + len(self.StateGraphs[k,j].nodes())) == len(self.WeightGraph.nodes())

            # Due to the fact that we're eliminating self/self transitions (to save computational time and memory), it's not unusual to have totally empty graphs if we run this
            # too early.  Don't panic; it just means we haven't established a proper state to state network yet.  It should line up nicely with what's in the assignment file, if
            # you check it.



            pi.new_operation('Determining correlation and eliminating correlated nodes...', len(start_pts))
            print("Okay but really?")
            for k in xrange(nstates):
                for j in xrange(nstates):
                    if k != j:
                    #for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
                        #iter_nodes = self.StateGraphs[k,j].nodes(data=True)
                        iter_nodes = self.WeightGraph.nodes(data=True)
                        niter_nodes = []
                        for i in iter_nodes:
                            if i[1]['iteration'] == 1:
                                niter_nodes.append(i)
                        # We have a list of the nodes that appear in the first iteration.  With this, we'll start going over the children, determining whether we have
                        # two or more children, then analysing their correlation length in the pcoord.  If there is a merge down the line, but the correlation length
                        # is longer than that, we'll trim all but one of the nodes (it doesn't particularly matter which one, really).
                        # Ideally, these correlated walkers should be switching in and out of the same bins, but it would be good to assert that.
                        # Actually, for an accurate count, we could probably just:
                        #   1. determine the true correlation length, which may be multiple tau, then remove nodes if they merge before the correlation length is 'over', or
                        #   2. determine if the correlation length is longer than 1 tau, then remove the node and add edges from its children to the remaining node.
                        # 2 has the distinct disadvantage of giving us less interesting information to work with in the future, however, although it is easier to implement.
                        # It would certainly be accurate in terms of the number of pathways, however, but that is no longer all we care about.
                        # Probably the best way to do this, in order to keep as much information as possible to play around with later, is to do the following:
                        # Create a correlation function.
                        # Run over the network, and every time there's a split, feed the resultant pcoord into the split (going into the pcoord of the children's children, if necessary)
                        # Record the resultant information in a list in a bin-tagged list, as well as on the node itself.
                        # Once we've run through the whole network, return to the beginning.
                        # Determine if the results of a split remerge before their correlation length has ended.
                        # If so, prune (or merge, essentially; add edges from removed nodes' children to remaining node).
                        # Generate the count matrix by using bin_assignments[0] and bin_assignments[-1] information and adding +1 to the respective information in the count matrix.
                        for i in niter_nodes:
                            #print(i)
                            #print(len(self.StateGraphs[k,j].neighbors(i[0])))
                            #print(self.StateGraphs[k,j].neighbors(i[0]))
                            #if len(self.StateGraphs[k,j].neighbors(i[0])) > 1:
                            if len(self.WeightGraph.neighbors(i[0])) > 1:
                                child_list = []
                                pcoord_child_list = nx.get_node_attributes(self.WeightGraph, 'pcoord')
                                #for children in self.StateGraphs[k,j].neighbors_iter(i[0]):
                                for children in self.WeightGraph.neighbors_iter(i[0]):
                                    # TEST
                                    #child_list.append(self.StateGraphs[k,j].node(children))
                                    child_list.append(children)
                                #print(pcoord_child_list[child_list[0]][:,0])
                                #print(pcoord_child_list[child_list[1]][:,0])
                                #print(child_list)
                                # Okay, so we have our pcoord...
                                #corr = np.array(signal.correlate(pcoord_child_list[child_list[0]][:,0], pcoord_child_list[child_list[1]][:,0], mode='full'))
                                # We can normalise our pcoord vectors before doing the cross correlation to get a normalised output; at 1, the vectors are correlated, and
                                # eventually, they'll be... well, not 1.
                                #a = np.fft.fft(pcoord_child_list[child_list[0]][:,0] / np.linalg.norm(pcoord_child_list[child_list[0]][:,0], axis=0))
                                #b = np.fft.fft(pcoord_child_list[child_list[1]][:,0] / np.linalg.norm(pcoord_child_list[child_list[1]][:,0], axis=0))
                                #a = pcoord_child_list[child_list[0]][:,0] / np.linalg.norm(pcoord_child_list[child_list[0]][:,0], axis=0)
                                #b = pcoord_child_list[child_list[1]][:,0] / np.linalg.norm(pcoord_child_list[child_list[1]][:,0], axis=0)
                                a = pcoord_child_list[child_list[0]][:,0]
                                b = pcoord_child_list[child_list[1]][:,0]
                                #print(a.sum())
                                print(a, b)
                                #corr = np.array(np.correlate(a, b, mode='full'))
                                coh = cohere(a, b, NFFT=5)[0]
                                print(coh)
                                plt.plot(coh)
                                #uconf = [1.96 / np.sqrt(len(corr))] * len(corr)
                                #lconf = [-1.96 / np.sqrt(len(corr))] * len(corr)
                                #plt.plot(uconf)
                                #plt.plot(lconf)
                                plt.show()
                                #x = range(0,len(pcoord_child_list[child_list[0]][:,0]))
                                #xcorr = np.arange(corr.size)
                                #lags = xcorr - (pcoord_child_list[child_list[0]][:,0].size-1)
                                #distancePerLag = (x[-1] - x[0])/float(len(x))
                                #print(len(corr))
                                #print(lags)
                                #print(distancePerLag)
                                #offsets = -lags*distancePerLag
                                #print(offsets)
                                # Okay, so it's unsurprising that the highest correlation value is set at time 0, I think, if I'm reading this correctly.  Not a big shock, given what we're dealing with.
                                # I should just see about converting it back into a function where it's indexed to the normal pcoord time steps, because I think it's symmetrical due to the nature of the 
                                # function and we can just normalise according to the uber peak.
                                # Wait!  No, I can't do that.  That would remove vital information.  However, I absolutely can normalise according to that peak, and then move away from that point and...
                                # wait... no.  Shouldn't be that symmetric?  Gotta think about this for a second.  Why are... oh, because I locked it to one piece of data in the 'lags' section, that's why.






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
