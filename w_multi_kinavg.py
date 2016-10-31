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
from westtools.dtypes import iter_block_ci_dtype as ci_dtype

#ci_dtype = numpy.dtype([('iter_start', n_iter_dtype),
#                        ('iter_stop', n_iter_dtype),
#                        ('expected', numpy.float64),
#                        ('ci_lbound', numpy.float64),
#                        ('ci_ubound', numpy.float64),
#                        ('corr_len', n_iter_dtype),
#                        ('variance', numpy.float64),
#                        ('stderrormean', numpy.float64)])

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
        # We no longer care about a lot of this.
        self.ntrials = 0
        self.nstates = 0
        self.kin_trial = {}
        self.west = {}
        self.niters = 0

    def add_args(self, parser):
        self.progress.add_args(parser)
        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('--output-file', default='m_kinavg.h5',
                            help='''The name of the output file to store results in.''')
        iogroup.add_argument('-k','--kinetics', default='kinavg.h5', 
                            help='''The name of the main .h5 file inside each simulation
                             directory''')


    def open_files(self):
        self.output_file = h5io.WESTPAH5File(self.output_file, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)

        opened_files = self.generate_file_list([self.west])
        self.kinH5 = opened_files[self.kin]
        # Just some temp things while I clean everything up...
        west_files = self.westH5
        # Determine max iteration ...

        # We can't really use the old method anymore, as we need to calculate rates in the bootstrap.
        # Ergo, we're going to load things like w_kinavg, but that's all.
        # We'll just load them up and store them internally, for the moment.

    def process_args(self, args):
        self.progress.process_args(args)
        self.output_file = args.output_file
        self.kin = args.kinetics
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
            # self.total_number_of_walkers()


            # Create a giant WEST.h5 file, separating the individual walkers, and renormalizing the weights.
            # It should then be compatible with existing toolsets.
            # Isn't really going to start with auxdata, but we'll add it in.

            #self.niters = 500
            pi.new_operation('Recreating...', self.niters)
            kinh5 = []
            for ifile, (key, west) in enumerate(self.westH5.iteritems()):
                kinh5.append(west)
                # self.niters = west.attrs['west_current_iteration'] - 1
            start_point = []

            # Initialize datasets.
            target_evol = numpy.zeros((len(start_pts), nstates), dtype=ci_dtype)
            flux_evol = numpy.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
            rate_evol = numpy.zeros((len(start_pts), nstates, nstates), dtype=ci_dtype)
            for iter in range(self.niters):
                ctarget = numpy.zeros((len(kinh5), nstates), dtype=ci_dtype)
                cflux = numpy.zeros((len(kinh5), nstates, nstates), dtype=ci_dtype)
                crate = numpy.zeros((len(kinh5), nstates, nstates), dtype=ci_dtype)
                # Load up, average, then call it a day.
                for ifile, kinetics in enumerate(kinh5):
                    ctarget[ifile, :] = kinetics['target_flux_evolution'][iter,:]
                    cflux[ifile, :, :] = kinetics['conditional_flux_evolution'][iter,:,:]
                    crate[ifile, :, :] = kinetics['rate_evolution'][iter,:,:]
                # First, we'll calculate the average...
                target_evol[iter,:]['expected'] = np.average(ctarget[:,:]['expected'], axis=0)
                flux_evol[iter,:,:]['expected'] = np.average(cflux[:,:,:]['expected'], axis=0)
                rate_evol[iter,:,:]['expected'] = np.average(crate[:,:,:]['expected'], axis=0)
                # ... then we'll sort the upper and lower CI.
                target_error = 1.96*np.std(ctarget[:,:]['expected'], axis=0)
                flux_error = 1.96*np.std(cflux[:,:,:]['expected'], axis=0)
                rate_error = 1.96*np.std(crate[:,:,:]['expected'], axis=0)
                target_evol[iter,:]['ci_lbound'] = target_evol[iter,:]['expected'] - target_error
                flux_evol[iter,:,:]['ci_lbound'] = flux_evol[iter,:,:]['expected'] - flux_error
                rate_evol[iter,:,:]['ci_lbound'] = rate_evol[iter,:,:]['expected'] - rate_error
                target_evol[iter,:]['ci_ubound'] = target_evol[iter,:]['expected'] + target_error
                flux_evol[iter,:,:]['ci_ubound'] = flux_evol[iter,:,:]['expected'] + flux_error
                rate_evol[iter,:,:]['ci_ubound'] = rate_evol[iter,:,:]['expected'] + rate_error
                # Then we'll get the start and stop iter.
                target_evol[iter,:]['iter_start'] = ctarget[:,:]['iter_start']
                flux_evol[iter,:,:]['iter_start'] = cflux[:,:,:]['iter_start']
                rate_evol[iter,:,:]['iter_start'] = crate[:,:,:]['iter_start']
                target_evol[iter,:]['iter_stop'] = ctarget[:,:]['iter_stop']
                flux_evol[iter,:,:]['iter_stop'] = cflux[:,:,:]['iter_stop']
                rate_evol[iter,:,:]['iter_stop'] = crate[:,:,:]['iter_stop']
                # For now, who cares about the others?
                # They are: corr_len, variance, stderrormean


                        

                pi.progress +=1 
        pi.new_operation('Writing to file...')
        df_ds = self.output_file.create_dataset('conditional_flux_evolution', data=flux_evol, shuffle=True, compression=9)
        tf_ds = self.output_file.create_dataset('target_flux_evolution', data=target_evol, shuffle=True, compression=9)
        #cp_ds = self.output_file.create_dataset('color_prob_evolution', data=pops.data, shuffle=True, compression=9)
        rate_ds = self.output_file.create_dataset('rate_evolution', data=rate_evol, shuffle=True, compression=9)
        
        for ds in (df_ds, tf_ds, rate_ds):
            self.stamp_mcbs_info(ds)
        ds_rate_evol = self.output_file.create_dataset('summary', data=summary, shuffle=True, compression = 9)

if __name__ == '__main__':
   WMultiKinetics().main()

