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
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

# directory locations are stored in a .yaml file with this format:
# ---
# PATHS: ['/path/to/simulation/1','/path/to/simulation/2',...,
# '/path/to/simulation/n']

class WMultiSimTool(WESTTool):
    prog ='w_multi_sim_tool'
    description = '''\
Calculate rate constant and associated errors by averaging over multiple 
simulations WE h5 files from each simulations are required to be stored 
in a .yaml format. (See "w_multli_sim_tool --help" for more information).
-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WESTTool,self).__init__()
        self.progress = ProgressIndicatorComponent()
        self.output_file = None
        self.input_file = None
        self.ntrials = 0
        self.nstates = 0
        self.niters = sys.maxint

    def add_args(self, parser):
        self.progress.add_args(parser)
        iogroup = parser.add_argument_group('input/output options')

        iogroup.add_argument('-o', '--output-file', default='multi_sim_output.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        iogroup.add_argument('-i','--input-file', default='simulation_directories.yaml', 
                            help='''File storing locations of .h5 files for each simulation in
                             yaml format (default: %(default)s).''')

    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_file, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)

    def process_args(self, args):
        self.progress.process_args(args)
        self.output_file = args.output_file
        self.input_file = args.input_file
        directory_dictionary = yaml.load(open(self.input_file))
        self.ntrials = len(directory_dictionary['PATHS'])
        for trial in range(self.ntrials):
            kinavg_file = h5py.File(directory_dictionary['PATHS'][trial] + '/kinavg.h5','r')
            self.nstates = len(kinavg_file['state_labels'])
            if len(kinavg_file['rate_evolution']) < self.niters:
                self.niters = len(kinavg_file['rate_evolution'])

    def go(self):
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Initializing')
            self.open_files()
            
            # We'll eventually store our data in the following structure. There is almost
            # certainly an opportunity to reduce memory use down the line, but we'll
            # worry about that once everything else is up and running.
            
            avgd_flux_evol = numpy.zeros((self.niters, self.nstates, self.nstates), dtype=ci_dtype)
            
            directory_list = yaml.load(open(self.input_file))['PATHS']
            
            # Data structure for storing average values is in an niter by nstates by nstates
            # matrix (essentially a loaf of bread, with each slice being an nstates by nstates
            # square, and the loaf contains niters slices). avg_loaf[100][0][1] will store
            # the average rate from state 0 to 1 at iteration 100.
            
            avg_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            
            # We'll use a similar structure to store the bounds on the confidence intervals
            lbound_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            ubound_loaf = numpy.zeros( (self.niters, self.nstates, self.nstates) )
            
            # Now, let's loop over each slice of bread.
            pi.new_operation('Calculating trial-average rates and confidence intervals...', self.niters)
            for iter in range(self.niters):
                # We're going to grab each slice of bread and cumulatively average the rates
                # between each state by just adding expected rate evolution data from each
                # trial, then normalizing in the end by the number of trials.
                for trial in range(self.ntrials):
                    kin_trial = h5py.File(directory_list[trial] + '/kinavg.h5')
                    avg_loaf[iter,:,:] = avg_loaf[iter,:,:] + kin_trial['rate_evolution']['expected'][iter,:,:]
                
                avg_loaf[iter,:,:] = numpy.divide(avg_loaf[iter,:,:],self.ntrials)
                sigma = numpy.zeros( (self.nstates, self.nstates) )

                # Now that we have the average values for the rates between states for
                # this iter, let's also go ahead and calculate the upper and lower bound
                # of the confidence interval.

                # Below we calculate the sample standard deviation among the trials:
                for trial in range(self.ntrials):
                    kin_trial = h5py.File(directory_list[trial] + '/kinavg.h5')
                    sigma = numpy.add( sigma, numpy.square( numpy.subtract( avg_loaf[iter,:,:],kin_trial['rate_evolution']['expected'][iter,:,:] ))/(self.ntrials-1))
                
                sigma = numpy.power(sigma,0.5)
                
                # We'll assume that we have enough independent trials to evoke the central
                # limit theorem, so we'll approximate the error using a Gaussian
                # distribution.
                
                for i in range(self.nstates):
                    for j in range(self.nstates):
                        if sigma[i,j] == 0:
                            lbound_loaf[iter,i,j] = 0
                            ubound_loaf[iter,i,j] = 0
                        else:
                            bounds = scipy.stats.norm.interval(0.95, loc = avg_loaf[iter,i,j], scale = sigma[i,j]/self.ntrials)
                            c_radius = (bounds[1]-bounds[0])/2
                            lbound_loaf[iter,i,j] = avg_loaf[iter,i,j] - c_radius
                            ubound_loaf[iter,i,j] = avg_loaf[iter,i,j] + c_radius
                
                        # Everything should be computed now so let's go ahead an wrap it all
                        # up in a neat little package.
                        avgd_flux_evol[iter]['expected'][i,j] = avg_loaf[iter,i,j]
                        avgd_flux_evol[iter]['ci_ubound'][i,j] = ubound_loaf[iter,i,j]
                        avgd_flux_evol[iter]['ci_lbound'][i,j] = lbound_loaf[iter,i,j]
                        
                
                pi.progress +=1 
            
            pi.new_operation('Writing to file...',1)
            ds_flux_evol = self.output_file.create_dataset('averaged_rate_evolution',data = avgd_flux_evol, shuffle = True, compression = 9)
            pi.progress +=1
            
if __name__ == '__main__':
   WMultiSimTool().main()
