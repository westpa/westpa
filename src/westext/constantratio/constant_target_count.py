# Modified from Josh's WExplore method
from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)

import numpy as np
import itertools

import westpa, west
from westpa import extloader
from westpa.yamlcfg import check_bool, ConfigItemMissing

class TargetRatio:
    def __init__(self, sim_manager, plugin_config):

        if not sim_manager.work_manager.is_master:
                return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        self.we_driver = sim_manager.we_driver
        self.priority = plugin_config.get('priority', 0)

        self.max_replicas = int(plugin_config.get('max_replicas'))
        self.states = plugin_config.get('state_definitions')
        self.state_to_trajectory = plugin_config.get('state_weights')
        if len(self.state_to_trajectory) == 1:
            self.state_to_trajectory = self.state_to_trajectory * len(self.states)

        # Register callback
        sim_manager.register_callback(sim_manager.prepare_iteration, self.post_we, self.priority)
        sim_manager.register_callback(sim_manager.pre_we, self.pre_we, self.priority)

    def pre_we(self):
        bin_mapper = self.system.bin_mapper
        segments = self.sim_manager.segments.values()

        final_pcoords = np.empty((len(segments), self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)

        for iseg, segment in enumerate(segments):
            final_pcoords[iseg] = segment.pcoord[-1,:]

        assignments = bin_mapper.assign(final_pcoords)

        state_bins = []
        #if self.states == 'None':
        if self.states == 'None':
            for i in (set(assignments)):
                state_bins.append(i)
        else:
            for i in self.states:
                state_bins.append(bin_mapper.assign([[i]]))

        # Now pull the 'ratio' (statically defined here, for the moment)
        active_bins = len(set(assignments))
        active_states = 0
        active_state_list = []
        for s_bin,t_bin in itertools.izip(state_bins, self.state_to_trajectory):
            if s_bin in assignments:
                active_bins += t_bin - 1
                active_states += 1
                active_state_list.append(s_bin)
        target_counts = np.empty((bin_mapper.nbins,), np.int_)
        bin_counts = int(np.floor(self.max_replicas / active_bins))
        if bin_counts < 1:
            westpa.rc.pstatus('')
            westpa.rc.pstatus('WARNING!:')
            westpa.rc.pstatus('Bins have been set to 0 walkers.  Enforcing a minimum of 1 walker per bin.  Consider adjusting your max walkers.  You WILL have more walkers this iteration.')
            westpa.rc.pstatus('')
            bin_counts = 1
        target_counts[...] = bin_counts

        extra = 0
        if active_states != 0:
            extra = np.floor((self.max_replicas - bin_counts*active_bins) / active_states)

        if extra < 0:
            extra = 0
        self.system.bin_target_counts = target_counts
        self.we_driver.bin_target_counts = target_counts
        for s_bin,t_bin in itertools.izip(state_bins, self.state_to_trajectory):
            self.system.bin_target_counts[s_bin] = (bin_counts * (t_bin)) + extra
            self.we_driver.bin_target_counts[s_bin] = (bin_counts * (t_bin)) + extra
        active_walkers = int(((active_bins)*bin_counts) + (extra*active_states))
        if active_walkers < self.max_replicas:
            self.system.bin_target_counts[active_state_list[np.random.randint(0,active_states)]] += (self.max_replicas - active_walkers)
            self.we_driver.bin_target_counts[active_state_list[np.random.randint(0,active_states)]] += (self.max_replicas - active_walkers)

        # Report stats
        westpa.rc.pstatus('-----------stats-for-next-iteration-')
        westpa.rc.pstatus('target counts: {}'.format(bin_counts))
        for ii_s_bin, s_bin in enumerate(state_bins):
            target = np.bincount(assignments)[s_bin]
            if target != 0 and self.states != 'None':
                #westpa.rc.pstatus('state {}, bin {} target counts: {}'.format(ii_s_bin, s_bin, np.bincount(assignments)[s_bin]))
                westpa.rc.pstatus('state {}, bin {} target counts: {}'.format(ii_s_bin, s_bin, self.we_driver.bin_target_counts[s_bin]))
        westpa.rc.pstatus('')
        westpa.rc.pflush()

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
                westpa.rc.pstatus('state {}, bin {} target counts: {}'.format(ii_s_bin, s_bin, np.bincount(assignments)[s_bin]))
        westpa.rc.pstatus('')
        westpa.rc.pflush()
