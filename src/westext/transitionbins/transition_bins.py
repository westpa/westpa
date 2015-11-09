# Modified from Josh's WExplore method
from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)

import numpy as np
import itertools

import westpa, west
from westpa import extloader
from westpa.yamlcfg import check_bool, ConfigItemMissing
from westpa.binning import RectilinearBinMapper

class TargetRatio:
    def __init__(self, sim_manager, plugin_config):

        if not sim_manager.work_manager.is_master:
                return

        # Define the bin types...
        function_mapper = {
                'RectilinearBinMapper': RectilinearBinMapper
                }
        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        self.we_driver = sim_manager.we_driver
        self.priority = plugin_config.get('priority', 0)
        n_iter = self.sim_manager.n_iter

        # This is non optional; do not disable.  There is no on-the-fly calculation possible.
        self.save_data = plugin_config.get('store_data', True)
        if self.sim_manager.n_iter == None or self.sim_manager.n_iter == 1:
            self.bin_boundaries = plugin_config.get('bin_boundaries', True)
            self.since_last_rebin = 0
            print('pulling from configuration...')
        else:
            self.bin_boundaries = self.data_manager.get_iter_group(n_iter)['transitionbins']['bin_boundaries']
            self.since_last_rebin = self.data_manager.get_iter_group(n_iter)['transitionbins']['since_last_rebin']
            print('pulling from .h5...')
        self.bin_mapper_type = function_mapper[plugin_config.get('bin_mapper', True)]
        self.bin_mapper = function_mapper[plugin_config.get('bin_mapper', True)]([self.bin_boundaries])
        self.bin_target_counts = plugin_config.get('target_counts', True)
        self.nbins = self.bin_mapper.nbins
        
        # Register callback
        #sim_manager.register_callback(sim_manager.prepare_iteration, self.post_we, self.priority)
        sim_manager.register_callback(sim_manager.pre_we, self.pre_we, self.priority)
        sim_manager.register_callback(sim_manager.prepare_iteration, self.prepare_new_iteration, self.priority)
        sim_manager.register_callback(sim_manager.post_we, self.prepare_new_iteration, self.priority)
        sim_manager.register_callback(sim_manager.prepare_run, self.prepare_new_iteration, self.priority)

    def normalize(self,m):
        # Taken shamelessly from Josh Adelman
        nm = m.copy()

        row_sum = m.sum(1)
        ii = np.nonzero(row_sum)[0]
        nm[ii,:] = m[ii,:] / row_sum[ii][:, np.newaxis]

        return nm

    def prepare_new_iteration(self):

        n_iter = self.sim_manager.n_iter
        bin_mapper = self.system.bin_mapper

        self.system.bin_mapper = self.bin_mapper_type([self.bin_boundaries])
        self.system.bin_target_counts = np.empty((bin_mapper.nbins,), np.int_)
        self.system.bin_target_counts[...] = 10
        print('preparing new iteration...')

    def pre_we(self):

        print('rebinning...')
        # This builds the transition matrix for the current set of bins.
        
        print(self.nbins)
        n_iter = self.sim_manager.n_iter
        bin_mapper = self.system.bin_mapper
        segments = self.sim_manager.segments.values()
        transition_matrix = np.zeros((bin_mapper.nbins, bin_mapper.nbins), dtype=np.float64)

        pcoords = np.empty((len(segments), 2, self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)
        assignments = np.empty((len(segments), 2), dtype=self.system.pcoord_dtype)

        self.system.bin_mapper = self.bin_mapper_type([self.bin_boundaries])
        self.system.bin_target_counts = np.empty((bin_mapper.nbins,), np.int_)
        self.system.bin_target_counts[...] = self.bin_target_counts

        for iseg, segment in enumerate(segments):
            pcoords[iseg,1] = segment.pcoord[-1,:]
            pcoords[iseg,0] = segment.pcoord[0,:]

        assignments[:,0] = bin_mapper.assign(pcoords[:,0])
        assignments[:,1] = bin_mapper.assign(pcoords[:,1])

        # Okay, so we've assigned the first and last time point of each walker into the assignment array.

        for iseg, segment in enumerate(segments):
            k = assignments[iseg,0]
            j = assignments[iseg,1]
            transition_matrix[k,j] += 1

        #transition_matrix = self.normalize(transition_matrix)

        # Write out that matrix!

        if self.save_data:
            # Create storage for ourselves
            with self.data_manager.lock:
                iter_group = self.data_manager.get_iter_group(n_iter)
                try:
                    del iter_group['transitionbins']
                except KeyError:
                    pass

        bin_iter_group = iter_group.create_group('transitionbins')
        bin_iter_group.create_dataset('transition_matrix', data=transition_matrix, compression=9)
        self.since_last_rebin += 1
        print(self.since_last_rebin)

        # Let's assume we're at the correct number of iterations to rebin.
        if self.since_last_rebin == 10:
            transition_matrix = np.empty((bin_mapper.nbins, bin_mapper.nbins), dtype=np.float64)
            low_value = max(1, n_iter-9)
            for iter in range(low_value,n_iter):
                transition_matrix = np.add(transition_matrix, self.data_manager.get_iter_group(iter)['transitionbins']['transition_matrix'][...])

            bin_object = {'old_index': 0, 'new_index': 0, 'active': True, 'modified': False}
            bins = {}
            boundary_width = np.diff(self.bin_boundaries[:-1])
            for ibin in range(0, bin_mapper.nbins-1):
                bins[ibin] = {'old_index': ibin, 'new_index': (ibin,0), 'active': True, 'modified': False, 'boundary_width': boundary_width[ibin], 'split': 1}
            transition_matrix = self.normalize(transition_matrix)
            # Let's do this in two rounds...
            index_adjustment = 0
            # Merge anything with no occupancy; pick the left or right bin at random
            for ibin in range(0, bin_mapper.nbins-1):
                if np.sum(transition_matrix[ibin,:]) == 0:
                    if np.sum(transition_matrix[:,ibin]) == 0:
                        random = np.random.randint(0,2)
                        if random == 0:
                            random = -1
                        if ibin == 0:
                            random = 1
                        if ibin == bin_mapper.nbins-2:
                            random = -1
                        bins[ibin+random]['boundary_width'] += bins[ibin]['boundary_width']
                        bins[ibin]['active'] = False
                        bins[ibin+random]['modified'] = False

            index_adjustment = 0
            for ibin in range(0, bin_mapper.nbins-1):
                if bins[ibin]['active'] == True:
                    transitions = np.where(transition_matrix[ibin,:] >= 0.15)[0]
                    # First, we check to see if any bins are too large...
                    if transitions.size != 0: 
                        loop_fired = True
                        largest_transition = np.sort(transitions)[-1]
                        if ibin == largest_transition and transitions.size != 1:
                            largest_transition = np.sort(transitions[-2])
                            if largest_transition == ibin:
                                break
                        if largest_transition > ibin and largest_transition != bin_mapper.nbins-1:
                            #if active_bins[ibin] == True:
                            if bins[ibin]['active'] == True:
                                if transition_matrix[ibin,largest_transition] >= 0.15:
                                    bins[ibin]['boundary_width'] += bins[ibin+1]['boundary_width']
                                    bins[ibin]['modified'] = True
                                    bins[ibin+1]['active'] = False
                        elif largest_transition < ibin:
                            if bins[ibin]['active'] == True:
                                if transition_matrix[ibin,largest_transition] >= 0.15:
                                    bins[ibin]['boundary_width'] += bins[ibin-1]['boundary_width']
                                    bins[ibin]['modified'] = True
                                    bins[ibin-1]['active'] = False

            # Adjust the indices...
            #index_adjustment = 0
            #for ibin in range(0, bin_mapper.nbins-1):
            #    bin = bins[ibin]
            #    if bin['active'] == True:
            #        bin['new_index'] -= index_adjustment
            #    else:
            #        index_adjustment += 1
            new_bins = {}
            #index_adjustment = 0
            for ibin in range(0, bin_mapper.nbins-1):
            #for ibin,bin in bins.iteritems():
                bin = bins[ibin]
                if bins[ibin]['active'] == True:
                    #bin['new_index'] += index_adjustment
                    if bins[ibin]['modified'] == False:
                        # Check to see if we're next to the 'sink' bin (that is, where things are stored so we don't oversample some boring region).
                        if transition_matrix[ibin,ibin] >= 0.90 and ibin != bin_mapper.nbins -1:
                            bound_width = bin['boundary_width']/2
                            bin['boundary_width'] = bound_width
                            bin['split'] = 3
                            index_adjustment += 2
                            new_bins[bin_mapper.nbins+index_adjustment-1] = {'old_index': [ibin,1], 'new_index': (ibin,1), 'active': True, 'modified': False, 'boundary_width': bound_width, 'split': 3}
                            new_bins[bin_mapper.nbins+index_adjustment] = {'old_index': [ibin,2], 'new_index': (ibin,2), 'active': True, 'modified': False, 'boundary_width': bound_width, 'split': 3}
                            # Now, modify the bins surrounding them...
                            for i,ad_bin in bins.iteritems():
                                if ad_bin['new_index'] == (bin['new_index'][0] - 1, ad_bin['split'] - 1):
                                    if ad_bin['active'] == True:
                                        print("Removing from previous bin!")
                                        ad_bin['boundary_width'] -= bound_width/2
                                # For those who haven't yet been modified by the index counter
                                if ad_bin['new_index'] == (bin['new_index'][0] + 1, 0):
                                    if ad_bin['active'] == True:
                                        print("Removing from next bin!")
                                        ad_bin['boundary_width'] -= bound_width/2
                        elif transition_matrix[ibin,ibin] >= 0.90 and ibin == bin_mapper.nbins - 1:
                            bound_width = bin['boundary_width']/3
                            bin['boundary_width'] = bound_width
                            bin['split'] = 2
                            index_adjustment += 1
                            new_bins[bin_mapper.nbins+index_adjustment] = {'old_index': (ibin,1), 'new_index': (ibin, 1), 'active': True, 'modified': False, 'boundary_width': bound_width*2, 'split': 2}
                        elif transition_matrix[ibin,ibin] >= 0.90 and ibin == 0:
                            bound_width = bin['boundary_width']/3
                            bin['boundary_width'] = bound_width*2
                            bin['split'] = 2
                            index_adjustment += 1
                            new_bins[bin_mapper.nbins+index_adjustment] = {'old_index': (ibin, 1), 'new_index': (ibin, 1), 'active': True, 'modified': False, 'boundary_width': bound_width, 'split': 2}
            active_bins = 0
            for ibin,bin in bins.iteritems():
                if bin['active'] == True:
                    active_bins += 1
            bin_widths = np.zeros(active_bins+len(new_bins.keys()))
            temp_list = {}
            print(bins)
            print(new_bins)
            for ibin,bin in bins.iteritems():
                if bin['active'] == True:
                    temp_list[bin['new_index']] = bin
            for ibin,bin in new_bins.iteritems():
                if bin['active'] == True:
                    temp_list[bin['new_index']] = bin
            current_bin = 0
            print(temp_list)
            for i in range(0, active_bins):
                print(i)
                try:
                    bin = temp_list[(i,0)]
                    split = bin['split']
                    for x in range(0, split):
                        if temp_list[(i,x)]['active'] == True:
                            bin_widths[current_bin] = temp_list[(i,x)]['boundary_width']
                            current_bin += 1
                            print(current_bin)
                except:
                    pass
                
            #for ibin,bin in bins.iteritems():
            #    if bin['active'] == True:
            #        bin_widths[bin['new_index']] = bin['boundary_width']
            #for ibin,bin in new_bins.iteritems():
            #        bin_widths[bin['new_index']] = bin['boundary_width']
            # A little funny business to handle the final index.
            print(bin_widths)
            self.bin_boundaries = [0] + list(np.cumsum(bin_widths)) + [self.bin_boundaries[-1]]
            self.since_last_rebin = 0
            print(self.bin_boundaries)

  
        bin_iter_group.create_dataset('since_last_rebin', data=(self.since_last_rebin))
        bin_iter_group.create_dataset('bin_boundaries', data=self.bin_boundaries, compression=9)
        self.system.bin_mapper = self.bin_mapper_type([self.bin_boundaries])
        self.system.bin_target_counts = np.empty((bin_mapper.nbins,), np.int_)
        self.system.bin_target_counts[...] = self.bin_target_counts
 
 
    def post_we(self):
        # There's a lot we don't care about here, such as coordinates, etc.  We're not changing the bin mapper, just the counts.
        segments = self.sim_manager.segments.values()
        bin_mapper = self.system.bin_mapper

        final_pcoords = np.empty((len(segments), self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)

        for iseg, segment in enumerate(segments):
            final_pcoords[iseg] = segment.pcoord[0,:]

        assignments = bin_mapper.assign(final_pcoords)

        # And ensure the ratio is set properly.  Let's define two particular bins as being 'state' bins, and adjust the count accordingly.
        # Actually, we'll need the assignments so that we know what bins are populated, now that I think about it.

        state_bins = []
        if self.states == 'None':
            for i in (set(assignments)):
                state_bins.append(i)
        else:
            for i in self.states:
                state_bins.append(bin_mapper.assign([[i]]))
        active_bins = len(set(assignments))
        active_states = 0
        for s_bin,t_bin in itertools.izip(state_bins, self.state_to_trajectory):
            if s_bin in assignments:
                active_bins += t_bin - 1
                active_states += 1
        target_counts = np.empty((bin_mapper.nbins,), np.int_)
        bin_counts = int(np.floor(self.max_replicas / active_bins))

        # Report level statistics
        westpa.rc.pstatus('')
        westpa.rc.pstatus('-----------stats-for-this-iteration-')
        westpa.rc.pstatus('target counts: {}'.format(bin_counts))
        #for ii_s_bin, s_bin in enumerate(state_bins):
        for ii_s_bin, s_bin in enumerate(state_bins):
            target = np.bincount(assignments)[s_bin]
            if target != 0 and self.states != 'None':
            #if target != 0:
                westpa.rc.pstatus('state {}, bin {} target counts: {}'.format(ii_s_bin, s_bin, np.bincount(assignments)[s_bin]))
        westpa.rc.pstatus('')
        westpa.rc.pflush()
