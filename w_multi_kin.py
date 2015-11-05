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
from matplotlib import pyplot

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

# directory locations are stored in a .yaml file with this format:
# ---
# PATHS: ['/path/to/simulation/1','/path/to/simulation/2',...,
# '/path/to/simulation/n']

class WMultiKinetics(WESTMultiTool):
    prog ='w_multi_kinetics'
    description = '''\
Calculate rate constant and associated errors by averaging over multiple 
simulations WE h5 files from each simulations are required to be stored 
in a .yaml format. (See "w_multi_sim_tool --help" for more information).
-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WESTTool,self).__init__()
        self.progress = ProgressIndicatorComponent()
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
        iogroup.add_argument('--non-markovian', action='store_true',
                            help='''Determine whether or not the output is from the non-Markovian
                            toolkit or not.''')
        iogroup.add_argument('-w','--west', default='west.h5', 
                            help='''The name of the main .h5 file inside each simulation
                             directory''')
        iogroup.add_argument('-a','--assign', default='assign.h5', 
                            help='''The name of the assignment .h5 file inside each simulation
                             directory''')
        iogroup.add_argument('-k','--kinetics', default='kintrace.h5', 
                            help='''The name of the kinetics .h5 file inside each simulation
                             directory.  Can be from non-Markovian or direct trace.''')
        iogroup.add_argument('--first_iter', default=1, 
                            help='''The iteration to start from.''')


    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_file, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        opened_files = self.generate_file_list([self.assign,self.west,self.kinetics])
        self.westH5 = opened_files[self.west]
        self.kineticsH5 = opened_files[self.kinetics]
        self.assignH5 = opened_files[self.assign]
        # Just some temp things while I cleane verything up...
        west_files = self.westH5
        directory_list = self.kineticsH5
        assign_list = self.assignH5
        # Determine max iteration ...
        iteration_list = []
        state_list = []
        for itrial, trial in enumerate(range(1, self.ntrials+1)):
            iteration_list.append(self.westH5[trial].attrs['west_current_iteration'] - 1)
            state_list.append(self.assignH5[trial].attrs['nstates'])
            if state_list[itrial] != state_list[itrial-1]:
                raise Exception('A different number of states have been used.')
        self.niters = min(iteration_list)
        self.nstates = state_list[0]

        # Create the datasets from scratch.  This avoids any messes with bootstrapping or older versions of the
        # non-Markovian toolkit.
        # We use this information to create a dataset object 'self.kin_trial', which is keyed to the trial ID
        for trial in range(1, self.ntrials+1):
            if self.non_markovian == False:
                h5 = directory_list[trial]
                assign = assign_list[trial]
                dataset = {}
                for key, value in h5.iteritems():
                    dataset[key] = value
                dataset['rate_evolution'] = numpy.zeros((self.niters, self.nstates, self.nstates), dtype=ci_dtype)
                dataset['conditional_flux_evolution'] = numpy.zeros((self.niters, self.nstates, self.nstates), dtype=ci_dtype)
                cfe = h5['conditional_fluxes']
                nstates = assign.attrs['nstates']
                nbins = assign.attrs['nbins']
                for istate in xrange(self.nstates):
                    for jstate in xrange(self.nstates):
                        flux = cfe[:,istate,jstate]
                        Pop = assign['labeled_populations'][:]
                        P = Pop[:,istate,:].sum(axis=1)
                        data = numpy.absolute(flux / P)
                        data[numpy.isnan(data)] = 0
                        data[numpy.isinf(data)] = 0
                        data_to_write = numpy.zeros(shape=self.niters)
                        cfe_to_write = numpy.zeros(shape=self.niters)
                        data = numpy.cumsum(data[self.first_iter:self.niters])
                        cfedata = numpy.cumsum(cfe[self.first_iter:self.niters,istate,jstate])
                        for ii, i in enumerate(data):
                            if ii != 0:
                                data[ii] = i / ii
                                cfedata[ii] = cfedata[ii] / ii
                        data_to_write[self.first_iter:self.niters] = data
                        cfe_to_write[self.first_iter:self.niters] = cfedata
                        dataset['rate_evolution'][:,istate,jstate]['expected'] = data_to_write
                        dataset['conditional_flux_evolution'][:,istate,jstate]['expected'] = cfe_to_write
                self.kin_trial[trial] = dataset
            else:
                h5 = h5py.File(directory_list[trial])
                dataset = {}
                # Should add all of our non-Markovian values in to the dictionary and make it seamless, regardless of which
                # toolkit we used.
                for key, value in h5.iteritems():
                    dataset[key] = value
                dataset['rate_evolution'] = numpy.zeros((self.niters, self.nstates, self.nstates), dtype=ci_dtype)
                cfe = h5['conditional_flux_evolution']
                for istate in xrange(self.nstates):
                    for jstate in xrange(self.nstates):
                        flux = cfe['expected'][:,istate,jstate]
                        P = h5['color_prob_evolution'][:,istate]
                        data = numpy.absolute(flux / P)
                        data[numpy.isnan(data)] = 0
                        data[numpy.isinf(data)] = 0
                        dataset['rate_evolution'][:,istate,jstate]['expected'] = data[:self.niters]
                self.kin_trial[trial] = dataset

    def monte_carlo(self, dataset, func=np.mean, alpha=.05, trials=1000):
        # These should eventually be replaced with Matt's Monte Carlo Bootstrap functions.
        trial_run = np.zeros(shape=trials)
        rand_array = dataset[np.random.randint(0, dataset.shape[0], size=(dataset.shape[0], trials))]
        trial_run[:] = func(rand_array, axis=0)
        trial_run = np.sort(trial_run)
        return_array = np.zeros((1), dtype=ci_dtype)
        return_array['expected'] = np.mean(trial_run)
        return_array['ci_lbound'] = trial_run[trials*(alpha/2)]
        return_array['ci_ubound'] = trial_run[trials*(1 - alpha/2)]
        return return_array

    def monte_carlo_pop(self, dataset, func=np.mean, alpha=.05, trials=1000):
        trial_run = []
        for i in xrange(0, trials):
            #for pull in dataset:
            trial_run.append(func(dataset[np.random.randint(0, dataset.shape[0], size=dataset.shape[0])]))
        trial_run = np.sort(trial_run)
        return_array = np.mean(trial_run)
        return return_array

    def process_args(self, args):
        self.progress.process_args(args)
        self.output_file = args.output_file
        self.non_markovian = args.non_markovian
        self.first_iter = args.first_iter - 1
        self.west = args.west
        self.assign = args.assign
        self.kinetics = args.kinetics
        self.sims = args.sims

    def total_number_of_walkers(self):
        self.total_walkers = [0]*self.niters
        for key,west in self.westH5.iteritems():
            # Sometimes, we're smaller or larger by one.  Hm.
            try:
                self.total_walkers[:] += west['summary'][:-1]['n_particles']
            except(ValueError):
                self.total_walkers[:] += west['summary'][:-1]['n_particles'][:len(self.total_walkers)]

    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            self.total_number_of_walkers()

            # We'll eventually store our data in the following structure. There is almost
            # certainly an opportunity to reduce memory use down the lines, but we'll
            # worry about that once everything else is up and running.

            avgd_rate_evol = numpy.zeros((self.niters, self.nstates, self.nstates), dtype=ci_dtype)
            avgd_flux_evol = numpy.zeros((self.niters, self.nstates, self.nstates), dtype=ci_dtype)
            avgd_color_prob = numpy.zeros((self.niters, self.nstates))
            avgd_state_prob = numpy.zeros((self.niters, self.nstates))

            avg_rate_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            avg_flux_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            avg_color_prob_loaf = numpy.zeros( (self.niters, self.nstates) )
            avg_state_prob_loaf = numpy.zeros( (self.niters, self.nstates) )

            # We'll use a similar structure to store the bounds on the confidence intervals
            lbound_rate_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            ubound_rate_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            lbound_flux_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            ubound_flux_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )

            # Now, let's loop over each slice of bread.
            pi.new_operation('Calculating trial-average rates and confidence intervals...', self.niters)
            for iter in range(self.niters):
                # We're going to grab each slice of bread and cumulatively average the rates
                # between each state by just adding expected rate evolution data from each
                # trial, then normalizing in the end by the number of trials.
                expected_rate = np.zeros((self.ntrials, self.nstates, self.nstates))
                expected_flux = np.zeros((self.ntrials, self.nstates, self.nstates))
                expected_state_prob = np.zeros((self.ntrials, self.nstates))
                expected_color_prob = np.zeros((self.ntrials, self.nstates))
                for i, (key,kin_trial) in enumerate(self.kin_trial.iteritems()):
                        expected_rate[i, :, :] = kin_trial['rate_evolution']['expected'][iter,:,:]
                        expected_flux[i, :, :] = kin_trial['conditional_flux_evolution']['expected'][iter,:,:]
                        try:
                            expected_state_prob[i,:,:] = kin_trial['state_prob_evolution'][iter,:]
                            expected_color_prob[i,:,:] = kin_trial['color_prob_evolution'][iter,:]
                        except:
                            pass
                for i in xrange(0, self.nstates):
                    for j in xrange(0, self.nstates):
                        if i != j:
                            avgd_rate_evol[iter,i,j] = self.monte_carlo(expected_rate[:,i,j], np.mean)
                            avgd_flux_evol[iter,i,j] = self.monte_carlo(expected_flux[:,i,j], np.mean)
                        avgd_color_prob[iter,i] = np.average(expected_color_prob[:,i])
                        avgd_state_prob[iter,i] = np.average(expected_state_prob[:,i])
                        

                pi.progress +=1 
        pi.new_operation('Writing to file...')
        ds_rate_evol = self.output_file.create_dataset('rate_evolution', data=avgd_rate_evol, shuffle=True, compression = 9)
        ds_rate_evol = self.output_file.create_dataset('conditional_flux_evolution', data=avgd_flux_evol, shuffle=True, compression = 9)
        ds_color_prob_evol = self.output_file.create_dataset('color_prob_evolution', data=avgd_color_prob, shuffle=True, compression = 9)
        ds_color_prob_evol = self.output_file.create_dataset('state_prob_evolution', data=avgd_state_prob, shuffle=True, compression = 9)
        ds_flux_evol = self.output_file.create_dataset('n_particles', data=self.total_walkers, shuffle=True, compression = 9)

if __name__ == '__main__':
   WMultiKinetics().main()

