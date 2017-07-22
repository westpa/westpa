# Modified from Josh's WExplore method
from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)

import numpy as np
import itertools

import westpa, west
from westpa import extloader
from westpa.yamlcfg import check_bool, ConfigItemMissing
from westpa.kinetics import RateAverager

class TargetRatio:
    def __init__(self, sim_manager, plugin_config):

        if not sim_manager.work_manager.is_master:
                return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        self.work_manager = sim_manager.work_manager
        self.we_driver = sim_manager.we_driver
        self.priority = plugin_config.get('priority', 0)
        self.N_sum = None

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

        self.automatic = True
        self.ratio_mode = False
        target_states = []
        for i in self.we_driver.target_states.values():
            target_states.append(bin_mapper.assign([i.pcoord])[0])

        if self.ratio_mode:
            state_bins = []
            #if self.states == 'None':
            if self.states == 'None':
                for i in (set(assignments)):
                    state_bins.append(i)
                self.state_to_trajectory = [1] * len(state_bins)
            else:
                for i in self.states:
                    state_bins.append(bin_mapper.assign([i]))

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
            active_walkers = int(((active_bins)*bin_counts))
            for s_bin,t_bin in itertools.izip(state_bins, self.state_to_trajectory):
                self.system.bin_target_counts[s_bin] = (bin_counts * (t_bin))
                self.we_driver.bin_target_counts[s_bin] = (bin_counts * (t_bin))
                #active_walkers += (bin_counts *(t_bin-1))
            #active_walkers = int(((active_bins)*bin_counts) + (active_states))
            if active_walkers < self.max_replicas:
                for i in xrange(0, self.max_replicas - active_walkers):
                    #rand = np.random.randint(0,active_states)
                    # Distribute amongst bins!  Why guarantee that we know best?
                    rand = np.random.randint(0,bin_mapper.nbins)
                    while rand not in assignments:
                        # Just pull a random number and make sure it's in the actual active bin list.
                        rand = np.random.randint(0,bin_mapper.nbins)
                    #self.system.bin_target_counts[active_state_list[rand]] += 1
                    #self.we_driver.bin_target_counts[active_state_list[rand]] += 1
                    self.system.bin_target_counts[rand] += 1
                    self.we_driver.bin_target_counts[rand] += 1
                    active_walkers += 1

        if self.automatic:
            # Create the transition rate averager, then, you know, use it.  Wait, does this do it properly?  This'll be a good time to check.
            n_iter = self.sim_manager.n_iter
            self.windowsize = 0.5
            eff_windowsize = int(n_iter * self.windowsize)
            averager = RateAverager(bin_mapper, self.system, self.data_manager, self.work_manager)
            averager.calculate(max(1, n_iter-eff_windowsize), n_iter+1, 1, 1)
            # Not currently sure if this will give us the rates directly without some manipulation.  However, we can just print it.
            #print(averager.average_rate.shape)
            #print(averager.average_rate)
            # Looks like those are indeed the rates.  Ergo, we want the self transition probability...
            # Set everything to flat sampling.
            active_bins = len(set(assignments)) - len(set(target_states))
            target_counts = np.empty((bin_mapper.nbins,), np.int_)
            bin_counts = int(np.floor(self.max_replicas / active_bins))
            target_counts[...] = bin_counts
            active_walkers = int(((active_bins)*bin_counts))
            accepted = 0
            if active_walkers < self.max_replicas:
                for i in xrange(0, self.max_replicas - active_walkers):
                    # Distribute amongst bins!  Why guarantee that we know best?
                    rand = np.random.randint(0,bin_mapper.nbins)
                    # Should REALLY properly handle recycling targets, now.
                    while rand not in assignments or rand in set(target_states):
                        # Just pull a random number and make sure it's in the actual active bin list.
                        rand = np.random.randint(0,bin_mapper.nbins)
                    target_counts[rand] += 1
                    active_walkers += 1
            
            # Now, we want to go through and adjust to ensure that our probabilities are okay...

            # Here's a function that will take a vector, p, and raise each element to the power of the corresponding element in vector N.
            # np.power(p, N)
            # First, check our distribution of transition probabilities and see how we're doing.  We want to minimize this.
            #p = averager.average_rate.diagonal().copy()
            # What if instead of minimizing the self transition probability, we MAXIMIZE the most unlikely transition probability?
            p = averager.average_rate.diagonal().copy()
            p = np.zeros(p.shape)
            p_i = np.zeros(p.shape)
            pcoord_var = np.zeros(p.shape)
            #print(averager.average_rate)
            for ii,i in enumerate(averager.average_rate):
                #print(i)
                try:
                    p_i[ii] = i[np.nonzero(i)].shape[0]
                    p[ii] = np.min(i[np.nonzero(i)])
                    # Assignments is the seg ID array corresponding to a bin assignment.
                    # Find the walkers that are in bin ii, then calculate their variance.
                    # We'll assume that our variance will increase until it's stable.
                    # Ergo, if the variance is larger than it was before, we'll assume the bin is 'done'.
                    pcoord_var[ii] = np.std(final_pcoords[np.where(assignments == ii)]) / np.average(final_pcoords[np.where(assignments == ii)])
                    #p[ii] = np.power(np.min(i[np.nonzero(i)]),p_i[ii])
                    # How many possible transitions are there?
                except:
                    pass
            p = 1 - p
            print(p)
            N = target_counts.copy()
            if self.N_sum == None:
                self.N_sum = np.zeros(N.shape[0])
                self.pcoord_var = np.zeros(p.shape)
            from scipy.stats import mode
            tp = np.power(p, (target_counts))
            tp_avg = np.average(tp)
            tp_err = np.std(tp) / tp_avg
            
            # Sort p such that it's in order of highest to lowest self transition probabilities...
            # First, sort lowest to highest, then reverse (the indexing).
            s = np.argsort(p)[::-1]
            #s = np.argsort(p)
            #s = np.arange(0, bin_mapper.nbins)

            # Iterate through the elements and see if we can't improve the current transition probabilities...
            print(s)
            for i in range(0, bin_mapper.nbins):
                # Check if the bin is active and not a recycling target.
                if s[i] in set(assignments) and s[i] not in set(target_states):
                    # Take from the bin with the LOWEST self transition probability...
                    #for j in reversed(range(0, bin_mapper.nbins)):
                    for j in reversed(range(0, bin_mapper.nbins)):
                        Continue = True
                        if s[j] in set(assignments):
                            while Continue:
                                # If we only have two transitions, that's likely the self transition and a 'backwards' transition.  Either way, probably
                                # want at least three valid transitions.  Enforcing this seems to cause problems, however, particularly with the edge region.
                                # Limit stealing from bin j if we have 10 times more sampling in the current bin i?
                                # Avoid divide by 0 errors.
                                # That part is actually unstable, though.
                                #if N[s[j]] > min(mode(p_i[np.nonzero(p_i)], axis=0)[0][0]*4, bin_counts) and float(self.N_sum[s[j]]) / max(float(np.sum(self.N_sum)),1) > .001:
                                # Let's make sure we give it a lot of chances to cross a barrier?
                                # When is the bin done sampling?  When the variance is relatively stable, maybe?
                                #if N[s[j]] > min(mode(p_i[np.nonzero(p_i)], axis=0)[0][0]*4, bin_counts) and self.N_sum[s[j]] > 100:
                                #if N[s[j]] > min(mode(p_i[np.nonzero(p_i)], axis=0)[0][0]*4, bin_counts) and pcoord_var[s[j]] < self.pcoord_var[s[j]]:
                                #if N[s[j]] > min(mode(p_i[np.nonzero(p_i)], axis=0)[0][0]*4, bin_counts):
                                #if N[s[j]] > 2:
                                if N[s[j]] > 2 and pcoord_var[s[j]] < self.pcoord_var[s[j]]:
                                    # We actually only care about the two bins, here...
                                    N[s[i]] += 1
                                    N[s[j]] -= 1
                                    # Check if we're good.  If yes, accept.  Otherwise, reject.
                                    new_tp = np.power(p, N)
                                    # We don't want to take too much from bins that have yet to be explored.  Let's try scoring whether or not we should actually take from bin s[j].
                                    #if np.std(tp[tp!=1])**2 - np.std(new_tp[new_tp!=1])**2 > 0 and (np.average(tp) - tp[s[i]]) > (np.average(new_tp) - tp[s[i]]):
                                    if np.std(tp[tp!=1])**2 - np.std(new_tp[new_tp!=1])**2 > 0:
                                        tp = new_tp
                                        accepted += 1
                                    else:
                                        N[s[i]] -= 1
                                        N[s[j]] += 1
                                        Continue = False

                                    #if N[s[j]] == 2:
                                    #if N[s[j]] == min(mode(p_i[np.nonzero(p_i)], axis=0)[0][0]*4, bin_counts):
                                    if N[s[j]] == 2:
                                        print("WHY AM I HERE")
                                        Continue = False
                                        break
                                else:
                                    Continue = False
                # One problem we have that is new bins are penalized, for some reason, with few walkers, likely due to their undersampled nature
                # (i.e., they have yet to sample rare transitions.  We should give them MORE walkers, not fewer).


                # To avoid this, we may wish to do a network deconstruction and see what bins are excluded from the strongly connected network
                # or possibly find bins for which there are fewer than average bin to bin transitions possible and give them a slight boost in 
                # the number of walkers.  That's not a bad idea, actually.  However, we'd need to be careful, as an edge bin may likely have fewer
                # possible pathways.

                else:
                    N[s[i]] = 0

                        
                        

            self.N_sum += N
            for ii,i in enumerate(self.pcoord_var):
                self.pcoord_var[ii] = max(i, pcoord_var[ii])
            # Boilerplate crap to make sure it doesn't fail out right now.
            self.system.bin_target_counts = N
            self.we_driver.bin_target_counts = N
            print(N)
            print(tp)
            print(self.pcoord_var)
            #print(np.power(p, N))
            bin_counts = N
            state_bins = []


        # Report stats
        westpa.rc.pstatus('-----------stats-for-next-iteration-')
        westpa.rc.pstatus('target counts: {}'.format(bin_counts))
        westpa.rc.pstatus('accepted rounds: {}'.format(accepted))
        for ii_s_bin, s_bin in enumerate(state_bins):
            target = np.bincount(assignments)[s_bin]
            if target != 0 and self.states != 'None':
            #if target != 0:
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
                state_bins.append(bin_mapper.assign([i]))
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
