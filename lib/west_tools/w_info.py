# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
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

from __future__ import print_function, division; __metaclass__ = type
import logging

# Let's suppress those numpy warnings.
import warnings
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=FutureWarning)

import sys, random, math
import numpy, h5py
from h5py import h5s

import westpa
from west.data_manager import weight_dtype, n_iter_dtype, seg_id_dtype
from westtools import (WESTMasterCommand, WESTTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent, BinMappingComponent)
from westpa import h5io

class WIWest(WESTSubcommand):
    subcommand = 'init'
    help_text = 'Pull information from a WESTPA HDF5 file without the configuration (rc) present.'
    description = '''\
'''
    def __init__(self, parent):
        super(WIWest,self).__init__(parent)
        
        # We're trying to initialize the west.h5 file, if available.
        # However, we can't guarantee that it exists.
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.binning = BinMappingComponent()
        self.data_manager = None
        
        self.output_filename = None
        # This is actually applicable to both.
        self.assignment_filename = None
        
        self.output_file = None
        self.assignments_file = None
        
        self.evolution_mode = None
        
        self.mcbs_alpha = None
        self.mcbs_acalpha = None
        self.mcbs_nsets = None

        # Now we're adding in things that come from the old w_kinetics
        self.do_compression = True
        
            
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        #subparsers = parser.add_subparsers(help='available commands')
        #info_parser = subparsers.add_parser('init', help='Display information about binning.')
        parser.add_argument('-n', '--n-iter', type=int, 
                                 help='''Consider initial points of segment N_ITER (default: current iteration).''')
        parser.add_argument('--detail', action='store_true',
                                 help='''Display detailed per-bin information in addition to summary
                                 information.''')
        suppress = ['--bins-from-system', '--bins-from-expr', '--bins-from-function', '--bins-from-file']
        self.binning.add_args(parser, suppress=suppress)
        #self.iter_range.include_args['iter_step'] = True
        #self.iter_range.add_args(parser)

        #iogroup = parser.add_argument_group('input/output options')
        #iogroup.add_argument('-a', '--assignments', default='assign.h5',
        #                    help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
        #                    (default: %(default)s).''')
        
        #iogroup.add_argument('-o', '--output', dest='output', default=self.default_output_file,
        #                    help='''Store results in OUTPUT (default: %(default)s).''')

    def process_args(self, args):
        # Open the data reader, which is necessary...
        self.data_reader.process_args(args)
        self.data_manager = self.data_reader.data_manager
        self.data_reader.open(mode='r')
        self.n_iter = getattr(args,'n_iter',None) or 1
        # We don't seem to know the bin hash, but heeeey.  If it the iter is 1, I think it's assumed/enforced
        # that the bin hash is what we initialized the simulation with.
        if self.n_iter == 1:
            self.binning.mapper_source_hash = self.data_manager.we_h5file['bin_topologies']['index']['hash'][0]
        self.binning.set_we_h5file_info(self.n_iter, self.data_reader)
        self.binning.process_args(args)
        #with self.data_reader:
        #    self.iter_range.process_args(args, default_iter_step=None)
        #if self.iter_range.iter_step is None:
            #use about 10 blocks by default
        #    self.iter_range.iter_step = max(1, (self.iter_range.iter_stop - self.iter_range.iter_start) // 10)
        
        #self.output_filename = args.output
        #self.assignments_filename = args.assignments

    def w_info(self):
        #print(dir(self.data_manager))
        #print(self.data_manager.we_h5file)
        #print(self.data_manager.current_iteration)
        #print(self.data_manager.get_iter_group(1).keys())
        #print(self.data_manager.get_iter_group(1).attrs.keys())
        #print(self.binning.mapper)
        # Okay, that seems to work, now.  We have the bin mapper, ergo we can...
        iter_group = self.data_manager.get_iter_group(self.n_iter)
        print(iter_group.keys())
        print(dir(self.binning.mapper))
        print(self.binning.mapper.boundaries)
        print(self.binning.mapper.labels)


    def go(self):
        self.w_info()

class WDirect(WESTMasterCommand, WESTTool):
    prog='w_direct'
    #subcommands = [AvgTraceSubcommand,AvgMatrixSubcommand]
    subcommands = [WIWest]
    subparsers_title = 'direct kinetics analysis schemes'

if __name__ == '__main__':
    WDirect().main()
