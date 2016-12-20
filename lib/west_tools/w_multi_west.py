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
            pi.new_operation('Recreating...', self.niters)
            westh5 = []
            for ifile, (key, west) in enumerate(self.westH5.iteritems()):
                westh5.append(west)
                self.niters = west.attrs['west_current_iteration'] - 1
            start_point = []
            for iter in range(self.niters):
                # We have the following datasets in each iteration:
                # ibstates, which aren't important.
                # pcoord
                # seg_index
                # wtgraph
                # wtgraph is going to be a little more complex to handle, but not too bad.
                iter += 1
                ifile = 0
                for west in westh5:
                    if iter == 1:
                        summary = west['summary'][...]
                    # We're going to (temporarily) assume we follow iter_00000 blah.
                    try:
                        seg_index = west['iterations/iter_{0:08d}'.format(iter)]['seg_index'][...]
                        pcoord = west['iterations/iter_{0:08d}'.format(iter)]['pcoord'][...]
                        wtgraph = west['iterations/iter_{0:08d}'.format(iter)]['wtgraph'][...]
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
                prev_start_point = start_point
                start_point = []
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
   WMultiKinetics().main()

