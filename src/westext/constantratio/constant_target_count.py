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

def flog(x, x0, a, b, c):
    y = a*np.log(b*(x-x0)) + c
    return y

# From https://gist.github.com/bistaumanga/6023692

'''Implementation and of K Means Clustering
Requires : python 2.7.x, Numpy 1.7.1+'''
#import numpy as np

def kMeans(X, K, maxIters = 10, plot_progress = None):

    centroids = X[np.random.choice(np.arange(len(X)), K), :]
    for i in range(maxIters):
        # Cluster Assignment step
        C = np.array([np.argmin([np.dot(x_i-y_k, x_i-y_k) for y_k in centroids]) for x_i in X])
        # Move centroids step
        centroids = [X[C == k].mean(axis = 0) for k in range(K)]
        if plot_progress != None: plot_progress(X, C, np.array(centroids))
    return np.array(centroids) , C

import sys
#plt.ion()

def show(X, C, centroids, keep = False):
    import time
    time.sleep(0.5)
    plt.cla()
    plt.plot(X[C == 0, 0], X[C == 0, 1], '*b',
         X[C == 1, 0], X[C == 1, 1], '*r',
         X[C == 2, 0], X[C == 2, 1], '*g')
    plt.plot(centroids[:,0],centroids[:,1],'*m',markersize=20)
    plt.draw()
    if keep :
        plt.ioff()
        plt.show()

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
        self.pcoord_var = None

        self.max_replicas = int(plugin_config.get('max_replicas'))
        self.states = plugin_config.get('state_definitions')
        self.state_to_trajectory = plugin_config.get('state_weights')
        # If we set it to automatic mode, then ignore any options relating to the ratio.
        #self.automatic = True if not plugin_config.get('automatic') else plugin_config.get('automatic')
        self.automatic = plugin_config.get('automatic')
        self.ratio_mode = False if self.automatic else True
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

        #self.automatic = True
        #self.ratio_mode = False
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
            accepted = 0

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
            # What if instead of minimizing the self transition probability, we MAXIMIZE the most unlikely transition probability?
            # We should get the shape, ultimately, but no one cares.
            p = averager.average_rate.diagonal().copy()
            p = np.zeros(p.shape)
            p_i = np.zeros(p.shape)
            index_p = np.zeros(p.shape)
            pcoord_var = np.zeros(p.shape)
            for ii,i in enumerate(averager.average_rate):
                #print(i)
                try:
                    p_i[ii] = i[np.nonzero(i)].shape[0]
                    p[ii] = np.min(i[np.nonzero(i)])
                    # This is just the bin ID that represents our lowest transition probability.
                    index_p[ii] = np.where(i == np.min(i[np.nonzero(i)]))[0][0]
                    # Assignments is the seg ID array corresponding to a bin assignment.
                    # Find the walkers that are in bin ii, then calculate their variance.
                    # We'll assume that our variance will increase until it's stable.
                    # Ergo, if the variance is larger than it was before, we'll assume the bin is 'done'.
                    pcoord_var[ii] = np.std(final_pcoords[np.where(assignments == ii)]) / np.average(final_pcoords[np.where(assignments == ii)])
                except:
                    pass
            # Transform from the probability of witnessing the rarest transition, into the probability of NOT witnessing the rarest transition.
            p = 1 - p
            print(p)
            N = target_counts.copy()
            try:
                I = self.I
            except:
                I = np.zeros((p.shape[0],2,self.max_replicas+1), dtype=float)
            if self.pcoord_var == None:
                try:
                    with self.data_manager.lock:
                        iter_group = self.data_manager.get_iter_group(n_iter-1)
                        # Pull in data from the west.h5 file.
                        self.pcoord_var = iter_group['ctc']['max_pcoord_var'][...]
                        # We don't actually seem to use this anymore, but we'll sort it later.
                        self.I = iter_group['ctc']['var_func_bin_count'][...]
                        I = self.I
                        self.N_sum = np.zeros(N.shape[0])
                        # Try pulling in the old bin counts..
                        # But, maybe we only want this is the number of occupied bins are the same, really.
                        #N = iter_group['ctc']['bin_target_counts'][...]
                        # It should probably work without doing that, actually.
                except:
                    self.N_sum = np.zeros(N.shape[0])
                    self.pcoord_var = np.zeros(p.shape)
            from scipy.stats import mode
            from scipy.optimize import curve_fit
            tp = np.power(p, (target_counts))
            tp_avg = np.average(tp)
            tp_err = np.std(tp) / tp_avg
            
            # Sort p such that it's in order of highest to lowest probability of not witnessing the rarest transitions.
            # First, sort lowest to highest, then reverse (the indexing).
            s = np.argsort(p)[::-1]
            #s = np.argsort(p)
            #s = np.arange(0, bin_mapper.nbins)

            # Iterate through the elements and see if we can't improve the current transition probabilities...
            print(s)
            # Default to two.
            try:
                assert minimum != None
            except:
                minimum = np.ones(p.shape) + 1
                #old_params = [(0,0]*p.shape[0]
            for i in range(0, bin_mapper.nbins):
                if i in set(assignments) and i not in set(target_states):
                    try:
                        #popt,pcov = curve_fit(flog, I[i,0][np.nonzero(I[i,0])],I[i,1][np.nonzero(I[i,0])], method='lm', maxfev=10000)
                        #x = np.linspace(0,self.max_replicas,self.max_replicas)
                        #y = flog(x, *popt)
                        X = np.vstack((I[i,0][np.nonzero(I[i,0])],I[i,1][np.nonzero(I[i,0])])).T
                        centroids, C = kMeans(X, K = 3, maxIters=10)
                        # Cut the minimum in half.  Is a threshold a thing that we want, here?
                        minimum[i] = max(int(np.floor(centroids[np.where(centroids[:,1] == np.max(centroids[:,1])),0]) / 4), 1)
                        if False:
                            if n_iter == 50:
                                from matplotlib import pyplot as plt
                                plt.scatter(I[i,0][np.nonzero(I[i,0])],I[i,1][np.nonzero(I[i,0])])
                                plt.plot(centroids[:,0],centroids[:,1],'*m',markersize=20)
                                show(X, C, centroids, True)
                                plt.show()
                    except Exception as e:
                        #print(e)
                        #minimum[i] = 2
                        pass
            print(minimum)
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
                                # Is the variance lower than the maximum?  Are we relatively close to the maximum variance?  If not, let that bin sample a bit more.
                                # What is the 'minimum?'  Calculate this based on a log fit.
                                #print(minimum[s[j]], N[s[j]])
                                # It's not so much a minimum as it is a 'don't take from me if I'm below this already'.
                                #if N[s[j]] > minimum[s[j]] and pcoord_var[s[j]] < self.pcoord_var[s[j]] and pcoord_var[s[j]] / self.pcoord_var[s[j]] > .5:
                                if N[s[j]] > minimum[s[j]] and pcoord_var[s[j]] < self.pcoord_var[s[j]]:
                                #if N[s[j]] > minimum[s[j]]:
                                    # We actually only care about the two bins, here...
                                    N[s[i]] += 1
                                    N[s[j]] -= 1
                                    # Check if we're good.  If yes, accept.  Otherwise, reject.
                                    new_tp = np.power(p, N)
                                    # We don't want to take too much from bins that have yet to be explored.  Let's try scoring whether or not we should actually take from bin s[j].
                                    #if np.std(tp[tp!=1])**2 - np.std(new_tp[new_tp!=1])**2 > 0 and (np.average(tp) - tp[s[i]]) > (np.average(new_tp) - tp[s[i]]):
                                    # SEEMS TO BE ALRIGHT!
                                    #if np.std(tp[tp!=1])**2 - np.std(new_tp[new_tp!=1])**2 > 0:
                                    # FAILURES
                                    #if np.std(np.log(tp[tp!=1])) - np.std(np.log(new_tp[new_tp!=1])) > 0:
                                    #if np.std(np.log(tp[tp!=1]))**2 - np.std(np.log(new_tp[new_tp!=1]))**2 > 0:

                                    # What if we just check that it minimizes the variance between the two bins?  That might work better.
                                    if np.std(tp[[s[j],s[i]]])**2 - np.std(new_tp[[s[j],s[i]]])**2 > 0:
                                        tp = new_tp
                                        accepted += 1
                                    else:
                                        N[s[i]] -= 1
                                        N[s[j]] += 1
                                        Continue = False

                                    #if N[s[j]] == 2:
                                    #if N[s[j]] == min(mode(p_i[np.nonzero(p_i)], axis=0)[0][0]*4, bin_counts):
                                    #if N[s[j]] == 2:
                                    #    print("WHY AM I HERE")
                                    #    Continue = False
                                    #    break
                                else:
                                    Continue = False
                else:
                    N[s[i]] = 0

                        
                        

            self.N_sum += N
            self.N = N
            for ii,i in enumerate(self.pcoord_var):
                self.pcoord_var[ii] = max(i, pcoord_var[ii])
                # Find if the element exists.  If so, replace it with the max.
                # First, see if it exists...
                #index = np.where(np.array(I[ii][0]) == N[ii])
                # We're doing this incorrectly, actually.
                #I[ii,1,int(N[ii])] = float(max(pcoord_var[ii], I[ii,1,int(N[ii])]))
                #I[ii,0,int(N[ii])] = int(N[ii])
                # Our pcoord variance is for the last iteration, not the current.
                I[ii,1,int(self.system.bin_target_counts[ii])] = float(max(pcoord_var[ii], I[ii,1,int(N[ii])]))
                I[ii,0,int(self.system.bin_target_counts[ii])] = int(self.system.bin_target_counts[ii])
            self.I = I
            # Boilerplate crap to make sure it doesn't fail out right now.
            self.system.bin_target_counts = N
            self.we_driver.bin_target_counts = N
            print(N)
            print(tp)
            print(self.pcoord_var)
            #print(np.power(p, N))
            bin_counts = N
            state_bins = []

            self.write_data = True
            if self.write_data:
                # Create storage for ourselves
                with self.data_manager.lock:
                    iter_group = self.data_manager.get_iter_group(n_iter)
                    try:
                        del iter_group['ctc']
                    except KeyError:
                        pass
        
                    ctc_iter_group = iter_group.create_group('ctc')
                    # This is the only one we currently care about from the averager.
                    ctc_iter_group.create_dataset('avg_rates', data=averager.average_rate, compression=4)
                    #wess_iter_group.create_dataset('unc_rates', data=averager.stderr_rate, compression=4)
                    # However, we DO care about our variance, our actual rate transition probabilities, and our target counts, etc.
                    ctc_iter_group.create_dataset('max_pcoord_var', data=self.pcoord_var, compression=4)
                    ctc_iter_group.create_dataset('pcoord_var', data=pcoord_var, compression=4)
                    ctc_iter_group.create_dataset('bin_target_counts', data=N, compression=4)
                    # The probably of witnessing an observation
                    ctc_iter_group.create_dataset('observ_prob', data=tp, compression=4)
                    # The rarest bin transition probabilities.
                    ctc_iter_group.create_dataset('rate_transition_prob', data=p, compression=4)
                    # We should probably save what the rarest transition actually IS, too.
                    ctc_iter_group.create_dataset('rarest_bin_transition_id', data=index_p, compression=4)
                    ctc_iter_group.create_dataset('number_of_possible_transitions', data=p_i, compression=4)
                    ctc_iter_group.create_dataset('var_func_bin_count', data=I, compression=4)
                    # We'll want to pull all this data, but that's for the next round of coding.


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
