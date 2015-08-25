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
import matplotlib
import argparse
from matplotlib import pyplot

from west.data_manager import weight_dtype, n_iter_dtype
from westtools import (WESTTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent)
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

class WMultiSimTool(WESTTool):
    prog ='w_multi_sim_tool'
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
        self.output_file = None
        self.input_file = None
        self.winput_file = None
        self.ntrials = 0
        self.nstates = 0
        self.kin_trial = {}
        self.west = {}
        self.niters = sys.maxint

    def add_args(self, parser):
        self.progress.add_args(parser)
        iogroup = parser.add_argument_group('input/output options')

        iogroup.add_argument('-o', '--output-file', default='multi_sim_output.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        iogroup.add_argument('-i','--input-file', default='simulation_directories.yaml', 
                            help='''File storing locations of .h5 files for each simulation in
                             yaml format (default: %(default)s).''')
        iogroup.add_argument('-wi','--winput-file', default='west_directories.yaml', 
                            help='''File storing locations of west.h5 files for each simulation in
                             yaml format (default: %(default)s).''')
        iogroup.add_argument('--non-markovian', action='store_true',
                            help='''Determine whether or not the output is from the non-Markovian
                            toolkit or not.''')
        # Ignore this.  It's garbage I put in here at 3 AM to help me effectively plot things for a meeting w/ Lillian.
        cogroup = parser.add_argument_group('temp options')
        cogroup.add_argument('--first-iter', type=int, default=1,
                            help='''Fraction of iterations to use in each window when running in ``cumulative`` mode.
                             The (1 - frac) fraction of iterations will be discarded from the start of each window.''')
        cogroup.add_argument('--block-size', type=int, default=1,
                            help='''Fraction of iterations to use in each window when running in ``cumulative`` mode.
                             The (1 - frac) fraction of iterations will be discarded from the start of each window.''')
        cogroup.add_argument('--tau', type=float, default=1.0,
                            help='''Fraction of iterations to use in each window when running in ``cumulative`` mode.
                             The (1 - frac) fraction of iterations will be discarded from the start of each window.''')



    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_file, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        west_files = yaml.load(open(self.winput_file))['PATHS']
        directory_list = yaml.load(open(self.input_file))['PATHS']
        for trial in range(self.ntrials):
            self.west[trial] = (h5py.File(west_files[trial]))
            if self.non_markovian == False:
                self.kin_trial[trial] = (h5py.File(directory_list[trial]))
            else:
                h5 = h5py.File(directory_list[trial])
                #flux = h5['conditional_flux_evolution']['expected']
                #P = h5['color_prob_evolution']
                dataset = {}
                # Should add all of our non-Markovian values in to the dictionary and make it seamless, regardless of which
                # toolkit we used.
                for key, value in h5.iteritems():
                    dataset[key] = value
                #print(dataset)
                dataset['rate_evolution'] = numpy.zeros((self.niters, self.nstates, self.nstates), dtype=ci_dtype)
                cfe = h5['conditional_flux_evolution']
                for istate in xrange(self.nstates):
                    for jstate in xrange(self.nstates):
                        flux = cfe['expected'][:,istate,jstate]
                        P = h5['color_prob_evolution'][:,istate]
                        #ii = numpy.where((P == 0.0))
                        data = numpy.absolute(flux / P)
                        data[numpy.isnan(data)] = 0
                        data[numpy.isinf(data)] = 0
                        #data[ii] = None
                        dataset['rate_evolution'][:,istate,jstate]['expected'] = data[:self.niters]
                self.kin_trial[trial] = dataset


    def process_args(self, args):
        self.progress.process_args(args)
        self.output_file = args.output_file
        self.input_file = args.input_file
        self.winput_file = args.winput_file
        self.non_markovian = args.non_markovian
        self.first_iter = args.first_iter
        self.block_size = args.block_size
        self.tau = args.tau
        directory_dictionary = yaml.load(open(self.input_file))
        self.ntrials = len(directory_dictionary['PATHS'])
        for trial in range(self.ntrials):
            kinavg_file = h5py.File(directory_dictionary['PATHS'][trial],'r')
            self.nstates = len(kinavg_file['state_labels'])
            if self.non_markovian == False:
                if len(kinavg_file['rate_evolution']) < self.niters:
                    self.niters = len(kinavg_file['rate_evolution'])
            else:
                if len(kinavg_file['conditional_flux_evolution']) < self.niters:
                    self.niters = len(kinavg_file['conditional_flux_evolution'])

    def total_number_of_walkers(self):
        self.total_walkers = [0]*self.niters
        for key,west in self.west.iteritems():
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

            #directory_list = yaml.load(open(self.input_file))['PATHS']

            # Data structure for storing average values is in an niter by nstates by nstates
            # matrix (essentially a loaf of bread, with each slice being an nstates by nstates
            # square, and the loaf contains niters slices). avg_loaf[100][0][1] will store
            # the average rate from state 0 to 1 at iteration 100.

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
                for key,kin_trial in self.kin_trial.iteritems():
                    #kin_trial = h5py.File(directory_list[trial])
                    #if self.non_markovian == False:
                        avg_rate_loaf[iter,:,:] += kin_trial['rate_evolution']['expected'][iter,:,:]
                        avg_flux_loaf[iter,:,:] += kin_trial['conditional_flux_evolution']['expected'][iter,:,:]
                        try:
                            avg_state_prob_loaf[iter,:] += kin_trial['state_prob_evolution'][iter,:]
                            avg_color_prob_loaf[iter,:] += kin_trial['color_prob_evolution'][iter,:]
                        except:
                            # Temp hack to get everything up and working.  Ultimately, we should just remove this if we don't have it.
                            pass
                    #else:
                    #    flux = kin_trial['conditional_flux_evolution']['expected']
                    #    P = kin_trial['color_prob_evolution']
                    #    dataset = {}
                    #    dataset['rate_evolution'] = numpy.zeros(flux.shape, dtype=ci_dtype)
                    #    for istate in xrange(self.nstates):
                    #        for jstate in xrange(self.nstates):
                    #            data = (flux[:,istate,jstate] / P[:, istate])
                    #            ii = numpy.where(P[:,istate] == 0.0)
                    #            data[ii] = None
                    #            dataset['rate_evolution'][:,istate,jstate] = data
                    #    kin_trial = dataset
                    #    avg_loaf[iter,:,:] = avg_loaf[iter,:,:] + kin_trial['rate_evolution']['expected'][iter,:,:]

                    

                avg_rate_loaf[iter,:,:] /= self.ntrials
                avg_flux_loaf[iter,:,:] /= self.ntrials
                avg_color_prob_loaf[iter,:] /= self.ntrials
                avg_state_prob_loaf[iter,:] /= self.ntrials
                sigma_rate = numpy.zeros( (self.nstates, self.nstates) )
                sigma_flux = numpy.zeros( (self.nstates, self.nstates) )

                # Now that we have the average values for the rates between states for
                # this iter, let's also go ahead and calculate the upper and lower bound
                # of the confidence interval.
                #for trial in range(self.ntrials):
                for key,kin_trial in self.kin_trial.iteritems():
                    #kin_trial = h5py.File(directory_list[trial])
                    #if self.non_markovian == True:
                    #    flux = kin_trial['conditional_flux_evolution']['expected']
                    #    P = kin_trial['color_prob_evolution']
                    #    dataset = {}
                    #    dataset['rate_evolution'] = numpy.zeros(flux.shape, dtype=ci_dtype)
                    #    for istate in xrange(self.nstates):
                    #        for jstate in xrange(self.nstates):
                    #            data = (flux[:,istate,jstate] / P[:, istate])
                    #            ii = numpy.where(P[:,istate] == 0.0)
                    #            data[ii] = None
                    #            dataset['rate_evolution'][:,istate,jstate] = data
                    #    kin_trial = dataset
                    #sigma_rate = numpy.add( sigma_rate, numpy.square( numpy.subtract( avg_rate_loaf[iter,:,:],kin_trial['rate_evolution']['expected'][iter,:,:] ))/(self.ntrials - 1))
                    #sigma_flux = numpy.add( sigma_flux, numpy.square( numpy.subtract( avg_flux_loaf[iter,:,:],kin_trial['conditional_flux_evolution']['expected'][iter,:,:] ))/(self.ntrials - 1))
                    sigma_rate = numpy.add( sigma_rate, numpy.square( numpy.subtract( avg_rate_loaf[iter,:,:],kin_trial['rate_evolution']['expected'][iter,:,:] )))
                    sigma_flux = numpy.add( sigma_flux, numpy.square( numpy.subtract( avg_flux_loaf[iter,:,:],kin_trial['conditional_flux_evolution']['expected'][iter,:,:] )))

                # Calculate the variance, assuming no standard error of the mean
                var_rate = sigma_rate / (self.ntrials-1)
                var_flux = sigma_flux / (self.ntrials-1)
                # Convert to standard error
                sigma_rate = numpy.power((sigma_rate/(self.ntrials-1)),0.5)
                sigma_flux = numpy.power((sigma_flux/(self.ntrials-1)),0.5)
                # We'll assume that we have enough independent trials to evoke the central
                # limit theorem, so we'll approximate the error using a Gaussian
                # distribution.

                for i in range(self.nstates):
                    for j in range(self.nstates):
                        if sigma_rate[i,j] == 0:
                            lbound_rate_loaf[iter,i,j] = 0
                            ubound_rate_loaf[iter,i,j] = 0
                        else:
                            bounds = scipy.stats.norm.interval(0.95, loc = avg_rate_loaf[iter,i,j], scale = sigma_rate[i,j] / numpy.sqrt(self.ntrials))
                            #bounds = scipy.stats.norm.interval(0.95, loc = avg_loaf[iter,i,j], scale = sigma[i,j])
                            c_radius = (bounds[1]-bounds[0])/2
                            lbound_rate_loaf[iter,i,j] = avg_rate_loaf[iter,i,j] - c_radius
                            ubound_rate_loaf[iter,i,j] = avg_rate_loaf[iter,i,j] + c_radius

                        if sigma_flux[i,j] == 0:
                            lbound_flux_loaf[iter,i,j] = 0
                            ubound_flux_loaf[iter,i,j] = 0
                        else:
                            bounds = scipy.stats.norm.interval(0.95, loc = avg_flux_loaf[iter,i,j], scale = sigma_flux[i,j] / numpy.sqrt(self.ntrials))
                            #bounds = scipy.stats.norm.interval(0.95, loc = avg_loaf[iter,i,j], scale = sigma[i,j])
                            c_radius = (bounds[1]-bounds[0])/2
                            lbound_flux_loaf[iter,i,j] = avg_flux_loaf[iter,i,j] - c_radius
                            ubound_flux_loaf[iter,i,j] = avg_flux_loaf[iter,i,j] + c_radius

                        # Everything should be computed now so let's go ahead an wrap it all
                        # up in a neat little package.
                        avgd_rate_evol[iter]['expected'][i,j] = avg_rate_loaf[iter,i,j] / self.tau
                        #avgd_rate_evol[iter]['variance'][i,j] = numpy.square(sigma_rate[i,j])
                        avgd_rate_evol[iter]['variance'][i,j] = var_rate[i,j]
                        avgd_rate_evol[iter]['stderrormean'][i,j] = sigma_rate[i,j] / numpy.sqrt(self.ntrials-1)
                        avgd_rate_evol[iter]['ci_ubound'][i,j] = ubound_rate_loaf[iter,i,j] / self.tau
                        avgd_rate_evol[iter]['ci_lbound'][i,j] = lbound_rate_loaf[iter,i,j] / self.tau
                        avgd_rate_evol[iter]['iter_start'][i,j] = self.first_iter
                        avgd_rate_evol[iter]['iter_stop'][i,j] = self.first_iter + (iter*self.block_size)

                        avgd_flux_evol[iter]['expected'][i,j] = avg_flux_loaf[iter,i,j] / self.tau
                        #avgd_flux_evol[iter]['variance'][i,j] = numpy.square(sigma_flux[i,j])
                        avgd_flux_evol[iter]['variance'][i,j] = var_flux[i,j]
                        avgd_flux_evol[iter]['stderrormean'][i,j] = sigma_flux[i,j] / numpy.sqrt(self.ntrials-1)
                        avgd_flux_evol[iter]['ci_ubound'][i,j] = ubound_flux_loaf[iter,i,j] / self.tau
                        avgd_flux_evol[iter]['ci_lbound'][i,j] = lbound_flux_loaf[iter,i,j] / self.tau
                        avgd_flux_evol[iter]['iter_start'][i,j] = self.first_iter
                        avgd_flux_evol[iter]['iter_stop'][i,j] = self.first_iter + (iter*self.block_size)
                
                    avgd_color_prob[iter,i] = avg_color_prob_loaf[iter,i]
                    avgd_state_prob[iter,i] = avg_state_prob_loaf[iter,i]

                pi.progress +=1 
            print(avgd_flux_evol)
        pi.new_operation('Writing to file...')
        ds_rate_evol = self.output_file.create_dataset('rate_evolution', data=avgd_rate_evol, shuffle=True, compression = 9)
        ds_rate_evol = self.output_file.create_dataset('conditional_flux_evolution', data=avgd_flux_evol, shuffle=True, compression = 9)
        ds_color_prob_evol = self.output_file.create_dataset('color_prob_evolution', data=avgd_color_prob, shuffle=True, compression = 9)
        ds_color_prob_evol = self.output_file.create_dataset('state_prob_evolution', data=avgd_state_prob, shuffle=True, compression = 9)
        ds_flux_evol = self.output_file.create_dataset('n_particles', data=self.total_walkers, shuffle=True, compression = 9)

if __name__ == '__main__':
   WMultiSimTool().main()

