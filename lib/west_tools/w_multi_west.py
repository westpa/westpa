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

import h5py
import numpy
import yaml
import sys
import scipy.stats
import logging
import re, os
import numpy, h5py
import numpy as np
import matplotlib
import argparse
import pickle
from matplotlib import pyplot
import logging
log = logging.getLogger(__name__)
from west import data_manager
from west.data_manager import weight_dtype, n_iter_dtype
from westtools import (WESTTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent, WESTMultiTool)
from westpa import h5io
#from westtools.dtypes import iter_block_ci_dtype as ci_dtype

ci_dtype = numpy.dtype([('iter_start', n_iter_dtype),
                        ('iter_stop', n_iter_dtype),
                        ('expected', numpy.float64),
                        ('ci_lbound', numpy.float64),
                        ('ci_ubound', numpy.float64),
                        ('corr_len', n_iter_dtype),
                        ('variance', numpy.float64),
                        ('stderrormean', numpy.float64)])

from random import randint

# directory locations are stored in a .yaml file with this format:
# ---
# PATHS: ['/path/to/simulation/1','/path/to/simulation/2',...,
# '/path/to/simulation/n']

# Straight up stolen from the data manager.  In the future, maybe I can just sort it by subbing in the appropriate values.
def get_bin_mapper(we_h5file,  hashval):
    '''Look up the given hash value in the binning table, unpickling and returning the corresponding
    bin mapper if available, or raising KeyError if not.'''

    # Convert to a hex digest if we need to
    try:
        hashval = hashval.hexdigest()
    except AttributeError:
        pass

    while True:
        # these will raise KeyError if the group doesn't exist, which also means
        # that bin data is not available, so no special treatment here
        try:
            binning_group = we_h5file['/bin_topologies']
            index = binning_group['index']
            pkl = binning_group['pickles']
        except KeyError:
            raise KeyError('hash {} not found. Could not retrieve binning group'.format(hashval))

        n_entries = len(index)
        if n_entries == 0:
            raise KeyError('hash {} not found. No entries in index'.format(hashval))

        chunksize = 1024

        for istart in xrange(0, n_entries, chunksize):
            chunk = index[istart:min(istart+chunksize, n_entries)]
            for i in xrange(len(chunk)):
                if chunk[i]['hash'] == hashval:
                    pkldat = bytes(pkl[istart+i, 0:chunk[i]['pickle_len']].data)
                    mapper = pickle.loads(pkldat)
                    log.debug('loaded {!r} from {!r}'.format(mapper, binning_group))
                    log.debug('hash value {!r}'.format(hashval))
                    return mapper

        raise KeyError('hash {} not found'.format(hashval))

class WMultiWest(WESTMultiTool):
    prog ='w_multi_west'
    description = '''\
Tool designed to combine multiple WESTPA simulations while accounting for
reweighting.  Test code thus far.
-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WESTTool,self).__init__()
        self.progress = ProgressIndicatorComponent()
        # We no longer care about a lot of this.
        self.ntrials = 0
        self.nstates = 0
        self.kin_trial = {}
        self.west = {}
        self.niters = 0

    def add_args(self, parser):
        self.progress.add_args(parser)
        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('--output-file', default='multi.h5',
                            help='''The name of the output file to store results in.''')
        iogroup.add_argument('-w','--west', default='west.h5', 
                            help='''The name of the main .h5 file inside each simulation
                             directory''')


    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_file, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)

        opened_files = self.generate_file_list([self.west])
        self.westH5 = opened_files[self.west]
        # Just some temp things while I clean everything up...
        west_files = self.westH5
        # Determine max iteration ...

        # We can't really use the old method anymore, as we need to calculate rates in the bootstrap.
        # Ergo, we're going to load things like w_kinavg, but that's all.
        # We'll just load them up and store them internally, for the moment.

    def process_args(self, args):
        self.progress.process_args(args)
        self.output_file = args.output_file
        self.west = args.west
        self.sims = args.sims

    def total_number_of_walkers(self):
        self.total_walkers = [0]*self.niters
        for key,west in self.westH5.iteritems():
            # Sometimes, we're smaller or larger by one.  Hm.
            try:
                self.total_walkers[:] += west['summary'][:-1]['n_particles']
            except(ValueError):
                self.total_walkers[:] += west['summary'][:-1]['n_particles'][:len(self.total_walkers)]

    class Segment():
        def __init__(self, weight=0, iteration=0, simid=0, recycled_in=0):
            self.weight = weight
            self.iteration = iteration
            self.simid = simid
            self.recycled_in = recycled_in

    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            self.total_number_of_walkers()

            # Create a giant WEST.h5 file, separating the individual walkers, and renormalizing the weights.
            # It should then be compatible with existing toolsets.
            # Isn't really going to start with auxdata, but we'll add it in.

            #self.niters = 500
            # Initialize data manager...
            self.data_manager = data_manager.WESTDataManager()
            pi.new_operation('Recreating...', self.niters)
            westh5 = []
            self.source_sinks = []
            self.n_sims = {}
            for ifile, (key, west) in enumerate(self.westH5.iteritems()):
                d = { 'west': west, 'wm': None, 'rt': None, 'remove_next_cycle': []}
                # We're getting the bin mapper, then setting the recycling target...
                binhash = west['iterations/iter_{0:08d}'.format(2)].attrs['binhash']
                bin_mapper = get_bin_mapper(west,binhash)
                print(bin_mapper.assign([[1.0]]))
                d['rt'] = bin_mapper.assign(west['tstates']['0']['pcoord'][...])[0]
                # We're going to make a set of source and sink states that we can iterate through, eventually.
                self.source_sinks.append(bin_mapper.assign(west['tstates']['0']['pcoord'][...])[0])
                # Keep a count of how many simulations for this particular recycling target we have...
                try:
                    self.n_sims[d['rt']] += 1
                except:
                    self.n_sims[d['rt']] = 1
                westh5.append(d)
                self.niters = west.attrs['west_current_iteration'] - 1
            start_point = []
            self.source_sinks = list(set(self.source_sinks))
            # We'll need a global list of walkers to add to and take care of during the next round of simulations, as well as the current one.
            # We'll organize it by source and sink states.
            self.past_iter = {}
            self.futr_iter = {}
            for i in self.source_sinks:
                self.past_iter[i] = []
                self.futr_iter[i] = []
            for iter in range(self.niters):
                # We have the following datasets in each iteration:
                # ibstates, which aren't important.
                # pcoord
                # seg_index
                # wtgraph
                # wtgraph is going to be a little more complex to handle, but not too bad.
                iter += 1
                ifile = 0
                # Determine how many simulations to append per west file.
                self.segments = {}
                for key,value in self.n_sims.iteritems():
                    self.segments[key] = int(np.floor(len(self.past_iter[key]) / value))

                for westdict in westh5:
                    # We have 2 key/value pairs in the dictionary: the opened HDF5 file, and WM, either None or a 'weight multiplier' matrix.
                    # Oh!  And a third one, now, too.  'rt', or the recycling target bin index.
                    west = westdict['west']
                    if iter == 1:
                        summary = west['summary'][...]
                    # We're going to (temporarily) assume we follow iter_00000 blah.
                    # We want and need to account for recycling, however.  That's going to require loading up the bin mapper and recycling index.
                    # We assume, for the moment, that the index doesn't change from iteration to iteration.  It'll be a one line change to make it happen,
                    # however.

                    # Try recreating the bin mapper...
                    binhash = west['iterations/iter_{0:08d}'.format(2)].attrs['binhash']
                    bin_mapper = get_bin_mapper(west,binhash)
                    # ... and done!  But a bit hackish, really.
                    # Okay, so we're going to make the assumption that the source of one simulation is the sink of the other.  That is, we're going to identify a set
                    # of recycling bins, then check whether any walkers have fallen into it.  If they have, we'll identify walkers in the OTHER simulation
                    # in the 'source' bin, and randomly add weight to them, then have it travel through by checking the parent weight multiplier.
                    # Copy in the previous weight matrix as the parent weight matrix for this simulation (west).
                    # Adjust the weights according to their multiplier.
                    # Identify recycling bins for this simulation.
                    # See if a walker has fallen into that state.
                    # If so, add the weight to a list to be randomly added to a new walker.
                    # During the next round, ensure the recycling bins do not equal each other.
                    # If they do not, take the weight and bin ID; find a walker within that bin, and add the weight appropriate, modified recycled weight.
                    # Calculate the new relative weight it should have, and modify its multiplier in the new weight matrix.
                    # We need to make sure we 'take out' the weight from the previous simulation too, though!  So, we need to keep a list of weight to (stochastically)
                    # remove from the current simulation during the next iteration, then update its relative weight matrix.

                    # So we need to define a segment object which we can push and pop onto a queue


                        # We have the bin hash and the recycling target.

                    # Check whether or not this exists... and try to change it, if it does.
                    #try:
                    #    if westdict['rt'] == None:
                    #        westdict['rt'] = bin_mapper.assign(w['tstates']['0']['pcoord'][0])
                    #except:
                        # Otherwise, this simulation doesn't have any recycling to worry about, so we're all good.
                    #    continue

                    try:
                        # Check whether the weight multiplication matrix exists for this code...
                        if westdict['wm'] == None:
                            seg_index = west['iterations/iter_{0:08d}'.format(iter)]['seg_index'][...]
                            westdict['wm'] = np.ones(seg_index.shape[0])
                            # ... and if not, create it.
                        else:
                            # Otherwise, transform the parent recycling matrix into this one.
                            # for the moment, just copy this code here (we'll rearrange later)
                            seg_index = west['iterations/iter_{0:08d}'.format(iter)]['seg_index'][...]
                            # Unfortunately, new segments have a negative index.  So we'll need to sort for those, first.
                            # We'll make a copy, then set all negatives to be 0.
                            # Afterwards, we'll return the index of where the parent_ids were negative,
                            # then set the multiplier to 1, as it should be.
                            unrecycled_parents = seg_index['parent_id'][...]
                            unrecycled_parents[np.where(unrecycled_parents < 0)[0]] = 0
                            westdict['wm'] = westdict['wm'][unrecycled_parents]
                            westdict['wm'][np.where(seg_index['parent_id'] < 0)[0]] = 1
                        seg_index = west['iterations/iter_{0:08d}'.format(iter)]['seg_index'][...]
                        pcoord = west['iterations/iter_{0:08d}'.format(iter)]['pcoord'][...]
                        wtgraph = west['iterations/iter_{0:08d}'.format(iter)]['wtgraph'][...]
                            # Let's reweight!
                            # First, we initialize to the correct weights based on prior iterations of recycling.
                        seg_index['weight'] *= westdict['wm']
                        # Now, we'll sort through any trajectories that need to be taken care of in terms of removing weight from this simulation...
                        # ... we need the source state, and the assignments.
                        source_state = self.source_sinks[np.where(np.array(self.source_sinks) != westdict['rt'])[0][0]]
                        assignments = bin_mapper.assign(pcoord[:,0,:])
                        for seg in westdict['remove_next_cycle']:
                            # We go through, pick a walker from the source state in the first time point, then remove the weight from it.
                            in_source = np.where(assignments == source_state)[0]
                            curr_weight = 0
                            max_weight = np.max(seg_index['weight'][in_source])
                            if max_weight > seg.weight:
                                while curr_weight <= seg.weight:
                                    remove_weight = in_source[randint(0,len(in_source)-1)]
                                    curr_weight = seg_index['weight'][remove_weight]
                                seg_index['weight'][remove_weight] -= seg.weight
                                # Now, we'll need to adjust the weight matrix entry for it.
                                # As we've already reweighted, it's easiest to just repull the damn thing and call it a day.
                                westdict['wm'][remove_weight] = seg_index['weight'][remove_weight] / west['iterations/iter_{0:08d}'.format(iter)]['seg_index']['weight'][remove_weight]
                            else:
                                # If there's no suitable walker, we just add it to everything there.
                                for remove_weight in in_source:
                                    seg_index['weight'][remove_weight] -= seg.weight / in_source.shape[0]
                                    # Now, we'll need to adjust the weight matrix entry for it.
                                    # As we've already reweighted, it's easiest to just repull the damn thing and call it a day.
                                    westdict['wm'][remove_weight] = seg_index['weight'][remove_weight] / west['iterations/iter_{0:08d}'.format(iter)]['seg_index']['weight'][remove_weight]


                        westdict['remove_next_cycle'] = []

                        # ... and then we'll handle adding in any new weight into our source state from other simulations.
                        #self.past_iter[i] = []
                        #self.futr_iter[i] = []
                        #self.n_sims[d['rt']] = 1
                        # We can pop items off the list; we just need to be sure of how many we're taking off.
                        # Okay, so we go through and pop off as many as we need to, assuming any exist.
                        # If we only have one left, just... add it.
                        if self.segments[source_state] != 0:
                            for i in range(0, self.segments[source_state]):
                                seg = self.past_iter[source_state].pop()
                                in_source = np.where(assignments == source_state)[0]
                                add_weight = in_source[randint(0,len(in_source)-1)]
                                seg_index['weight'][add_weight] += seg.weight
                                # Now, we'll need to adjust the weight matrix entry for it.
                                # As we've already reweighted, it's easiest to just repull the damn thing and call it a day.
                                westdict['wm'][add_weight] = seg_index['weight'][add_weight] / west['iterations/iter_{0:08d}'.format(iter)]['seg_index']['weight'][add_weight]
                        elif self.segments[source_state] == 0 and len(self.past_iter[source_state]) > 0:
                            for i in range(0, len(self.past_iter[source_state])):
                                seg = self.past_iter[source_state].pop()
                                in_source = np.where(assignments == source_state)[0]
                                add_weight = in_source[randint(0,len(in_source)-1)]
                                seg_index['weight'][add_weight] += seg.weight
                                # Now, we'll need to adjust the weight matrix entry for it.
                                # As we've already reweighted, it's easiest to just repull the damn thing and call it a day.
                                westdict['wm'][add_weight] = seg_index['weight'][add_weight] / west['iterations/iter_{0:08d}'.format(iter)]['seg_index']['weight'][add_weight]
                        if self.segments[source_state] > 1 and len(self.past_iter[source_state]) > 0:
                            seg = self.past_iter[source_state].pop()
                            in_source = np.where(assignments == source_state)[0]
                            add_weight = in_source[randint(0,len(in_source)-1)]
                            seg_index['weight'][add_weight] += seg.weight
                            # Now, we'll need to adjust the weight matrix entry for it.
                            # As we've already reweighted, it's easiest to just repull the damn thing and call it a day.
                            westdict['wm'][add_weight] = seg_index['weight'][add_weight] / west['iterations/iter_{0:08d}'.format(iter)]['seg_index']['weight'][add_weight]



                        # Then we'll handle our own recycling events, and add them to the list of things to handle in the next simulation.
                        assignments = bin_mapper.assign(pcoord[:,-1,:])
                        try:
                            in_sink = np.where(assignments[:] == westdict['rt'])[0]
                        except:
                            in_sink = []
                        for iseg in in_sink:
                            # We just add it to the list of stuff to remove next iteration.
                            westdict['remove_next_cycle'].append(self.Segment(weight=seg_index['weight'][iseg]))
                            self.futr_iter[westdict['rt']].append(self.Segment(weight=seg_index['weight'][iseg]))
                        

                        if ifile == 0:
                            mseg = seg_index
                            mpco = pcoord
                            mwtg = wtgraph

                            start_point.append(0)
                        if ifile != 0:
                            #print(mseg.shape, seg_index.shape, ifile)
                            #print(mpco.shape, pcoord.shape, ifile)
                            #print(mwtg.shape, wtgraph.shape, ifile)
                            if iter != 1:
                                addition = prev_start_point[ifile]
                            else:
                                addition = mseg.shape[0]
                            seg_index['parent_id'][np.where(seg_index['parent_id'] >= 0)] += addition
                            seg_index['parent_id'][np.where(seg_index['parent_id'] < 0)] -= addition
                            seg_index['wtg_offset'] += mwtg.shape[0]
                            start_point.append(mseg.shape[0])
                            wtgraph += mwtg.shape[0]
                            mseg = np.concatenate((mseg, seg_index))
                            mpco = np.concatenate((mpco, pcoord))
                            mwtg = np.concatenate((mwtg, wtgraph))
                        ifile += 1
                    except:
                        continue
                # Make a real copy to use in the next iteration.
                self.past_iter = self.futr_iter.copy()
                prev_start_point = start_point
                start_point = []
                # This is... maybe wrong, actually?  Or at least, it's not ALL that is required for normalizing things.
                # We need to weight everything by 1/N, then just normalize if that normalization was wrong.  Keep the relative weights sane.
                # ... or actually, no, that's fine, nevermind, what's wrong with me?  But we'll leave it in for now.
                mseg['weight'] /= mseg['weight'].sum()
                summary['n_particles'][iter-1] = mseg.shape[0]
                summary['norm'][iter-1] = mseg['weight'].sum()
                summary['min_seg_prob'][iter-1] = min(mseg['weight'])
                summary['max_seg_prob'][iter-1] = max(mseg['weight'])

                curr_iter = self.output_file.create_group('iterations/iter_{0:08d}'.format(iter))
                curr_iter.attrs['n_iter'] = iter
                ds_rate_evol = curr_iter.create_dataset('wtgraph', data=mwtg, shuffle=True, compression = 9)
                ds_rate_evol = curr_iter.create_dataset('seg_index', data=mseg, shuffle=True, compression = 9)
                ds_rate_evol = curr_iter.create_dataset('pcoord', data=mpco, shuffle=True, compression = 9)
                del mseg, mpco, mwtg

                        

                pi.progress +=1 
        pi.new_operation('Writing to file...')
        ds_rate_evol = self.output_file.create_dataset('summary', data=summary, shuffle=True, compression = 9)
        self.output_file.attrs['west_current_iteration'] = self.niters
        self.output_file.attrs['west_file_format_version'] = 7
        self.output_file.attrs['west_iter_prec'] = 8
        self.output_file.attrs['westpa_fileformat_version'] = 7
        self.output_file.attrs['westpa_iter_prec'] = 8

if __name__ == '__main__':
   WMultiWest().main()

