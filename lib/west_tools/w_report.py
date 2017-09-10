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
import itertools

# Let's suppress those numpy warnings.
import warnings
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=FutureWarning)

import sys, random, math
import numpy, h5py
import numpy as np
from h5py import h5s
import types

import westpa
from west.data_manager import weight_dtype, n_iter_dtype, seg_id_dtype
from westtools import (WESTMasterCommand, WESTTool, WESTDataReader, IterRangeSelection, WESTSubcommand,
                       ProgressIndicatorComponent, BinMappingComponent)
from westpa import h5io

class WIReport(WESTSubcommand):
    #subcommand = 'init'
    help_text = 'Pull information from a WESTPA HDF5 file with or without the configuration (rc) present.'
    description = '''\
'''
    def __init__(self, parent):
        super(WIReport,self).__init__(parent)
        
        # We're trying to initialize the west.h5 file, if available.
        # However, we can't guarantee that it exists.
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.binning = BinMappingComponent()
        self.data_manager = None
            
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        #subparsers = parser.add_subparsers(help='available commands')
        #info_parser = subparsers.add_parser('init', help='Display information about binning.')
        parser.add_argument('-n', '--n-iter', type=int, 
                                 help='''Consider initial points of segment N_ITER (default: current iteration).''')
        parser.add_argument('--detail', action='store_true',
                                 help='''Display detailed per-bin information in addition to summary
                                 information.''')
        parser.add_argument('-d', '--data', type=str, default=None,
                                 help='''The list of data to output.''')
        parser.add_argument('-l', '--line', action='store_true',
                                 help='''Report as single line.  Otherwise, use multiline key: value report.''')
        parser.add_argument('-s', '--separator', type=str,
                                 help='''Separator (delimiter) for use in line-mode.''')
        suppress = ['--bins-from-system', '--bins-from-expr', '--bins-from-function', '--bins-from-file']
        self.binning.add_args(parser, suppress=suppress)
        #self.iter_range.include_args['iter_step'] = True
        #self.iter_range.add_args(parser)

    def process_args(self, args):
        # Open the data reader, which is necessary...
        self.data_reader.process_args(args)
        self.data_manager = self.data_reader.data_manager
        self.data_reader.open(mode='r')
        self.n_iter = getattr(args,'n_iter', None) or 1
        self.data = getattr(args, 'data', None)
        if self.data is not None:
            self.data = self.data.split(' ')
        # We don't seem to know the bin hash, but heeeey.  If it the iter is 1, I think it's assumed/enforced
        # that the bin hash is what we initialized the simulation with.
        if self.n_iter == 1:
            self.binning.mapper_source_hash = self.data_manager.we_h5file['bin_topologies']['index']['hash'][0]
        self.binning.set_we_h5file_info(self.n_iter, self.data_reader)
        self.binning.process_args(args)
        self.line = getattr(args, 'line')
        self.args = args
        self.separator = getattr(args, 'separator')

    def report_single(self):
        west = self.WESTInfoBlob(self.data_reader, self.binning, self.args, self.n_iter)
        self.print_report(west, line=self.line, args=self.data, s2=self.separator)

    def report_range(self):
        west = self.WESTInfoBlob(self.data_reader, self.binning, self.args, self.n_iter)
        if self.separator == None:
            self.separator = ';'
        for i in range(1, self.n_iter+1):
            west.n_iter = i
            self.print_report(west, line=True, args=self.data, s2=self.separator)

    def lprint(self, k, l, ls, s=':'):
        import types, string
        if type(l) == types.ListType or type(l) == types.TupleType:
            #print(str(k).rjust(ls), s)
            print(string.capwords(str(k).replace('_',' ').replace('-',' ').replace('.', ' ')).rjust(ls), s)
            for i in l:
                self.lprint(k=' ',l=i,ls=ls,s=' ')
        elif type(l) == types.DictType:
            print(string.capwords(str(k).replace('_',' ').replace('-',' ').replace('.', ' ')).rjust(ls), s)
            label_size = 0
            for key, value in l.iteritems():
                label_size = max(label_size, len(key))
            for key, i in l.iteritems():
                self.lprint(k=key,l=i,ls=ls+label_size+len(s),s=':')
        else:
            print(string.capwords(str(k).replace('_',' ').replace('-',' ').replace('.', ' ')).rjust(ls), s, str(l).ljust(20))

    def deepgetattr(self, obj, attr):
        # We check to see if the object is callable.  If so, do it.  Why not?  Otherwise, no.
        # First, remove any function calls...
        try:
            # Not the first time we're doing this in WESTPA code.
            return eval(attr, {'__builtins__': {}}, obj.__dict__)
        except:
            return reduce (getattr, attr.split('.'), obj)

    def print_report(self, west, line=False, s1=':', s2=';', args=None):
        # Here, we're just going to print out these quantities...
        # We'll want to put in some more appropriate formatting eventually, but
        if args == None:
            args = ['n_iter', 'aggregate_walkers', 'recycling_events', 'tstates', 'mapper.labels', 'mapper.boundaries']
        output = {}
        label_size = 0
        for arg in args:
            try:
                if type(self.deepgetattr(west,arg)) == types.GeneratorType:
                    output[arg] = []
                    for t in self.deepgetattr(west,arg):
                        output[arg].append(t)
                else:
                    output[arg] = self.deepgetattr(west,arg)
                label_size = max(label_size, len(arg))
            except Exception as e:
                # Not an attribute.  Ergo, just... throw in an error, and go from there.
                output[arg] = 'ERROR - ' + str(e)
        if line == False:
            for arg in args:
                k = arg
                v = output[k]
                self.lprint(k,v,label_size, s1)
        else:
            print(s2.join(str(output[arg]) for arg in args))

    class WESTInfoBlob():
        '''
        A little class to ease pulling data in from WESTPA in a 'reporter friendly' way.
        Nothing too large, mind you.  But it should utilize the data manager and h5file without the need
        for the rc, and it should be able to return data that can be utilized in reporting.
        '''
        def __init__(self, data_reader, binning, args, n_iter = 1):
            # west is the west h5file from the data reader.
            # We should be able to change the iteration, if desired.
            self.data_reader = data_reader
            self.data_manager = self.data_reader.data_manager
            self.we_h5file = self.data_manager.we_h5file
            self.binning = binning
            self.mapper = None
            self.args = args
            self.iter_group = None
            self.n_iter = n_iter

        @property
        def n_iter(self):
            return self._n_iter

        @n_iter.setter
        def n_iter(self, n_iter):
            self._n_iter = n_iter
            if self._n_iter == 1:
                self.binning.mapper_source_hash = self.data_manager.we_h5file['bin_topologies']['index']['hash'][0]
            self.binning.set_we_h5file_info(self._n_iter, self.data_reader)
            self.binning.process_args(self.args)
            self.iter_group = self.data_manager.get_iter_group(self._n_iter)
            self.mapper = self.binning.mapper

        @property
        def tstates(self):
            if 'tstates' in self.iter_group.keys():
                for label,pcoord in itertools.izip(self.iter_group['tstates']['index'],self.iter_group['tstates']['pcoord']):
                    # We want the bin id, as well.  Call the bin mapper assign function.
                    binid = self.mapper.assign([pcoord])
                    yield {'Label': label,'Progress Coordinate': pcoord, 'Bin ID': binid}

        @property
        def recycling_events(self):
            # Just go through and sum up the aggregate number of events.
            self._recycling_events = 0
            for i in range(1, self.n_iter+1):
                iter_group = self.data_manager.get_iter_group(i)
                self._recycling_events += len(np.where(iter_group['seg_index']['endpoint_type'] == 3)[0])
            return self._recycling_events

        @property
        def aggregate_walkers(self):
            self._aggregate_walkers = self.data_manager.we_h5file['summary']['n_particles'][:self.n_iter].sum()
            return self._aggregate_walkers


class WISingle(WIReport):
    subcommand = 'single'
    help_text = 'Pull information from a WESTPA HDF5 file with or without the configuration (rc) present.'
    description = '''\
'''
    def go(self):
        self.report_single()

class WIRange(WIReport):
    subcommand = 'range'
    help_text = 'Pull information from a WESTPA HDF5 file with or without the configuration (rc) present.'
    description = '''\
'''
    def go(self):
        self.report_range()

class WDirect(WESTMasterCommand, WESTTool):
    prog='w_direct'
    #subcommands = [AvgTraceSubcommand,AvgMatrixSubcommand]
    subcommands = [WISingle, WIRange]
    subparsers_title = 'direct kinetics analysis schemes'

if __name__ == '__main__':
    WDirect().main()
