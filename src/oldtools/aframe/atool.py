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



import logging

log = logging.getLogger(__name__)

import os, sys
import h5py
import west, oldtools

from oldtools.aframe import AnalysisMixin, ArgumentError

class WESTAnalysisTool:
    def __init__(self):
        super(WESTAnalysisTool,self).__init__()
        # Whether a west.cfg is required to run a program based on this tool
        self.config_required = False
        
        # Analysis HDF5 filename and object
        self.anal_h5name = None
        self.anal_h5file = None
        
        # Whether this is being used in a brute analysis
        self.bf_mode = False
        
        # A way to override some arguments on a per-mixin basis without having to subclass
        # (messy, but it doesn't seem crucial enough so far to make it cleaner)
        self.include_args = {} 

    def add_args(self, parser, upcall = True):
        '''Add arguments to a parser common to all analyses of this type.'''
        if upcall:
            try:
                upfunc = super(WESTAnalysisTool,self).add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
                
        group = parser.add_argument_group('general analysis options')
        group.add_argument('-A', '--analysis-file', dest='anal_h5name', metavar='H5FILE', default='analysis.h5',
                            help='Store intermediate and final results in H5FILE (default: %(default)s).')
        
    def process_args(self, args, upcall = True):
        self.anal_h5name = args.anal_h5name
        
        if upcall:
            try:
                upfunc = super(WESTAnalysisTool,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)        

    def open_analysis_backing(self):
        if self.anal_h5file is None:
            self.anal_h5file = h5py.File(self.anal_h5name)
    
    def close_analysis_backing(self):
        try:
            self.anal_h5file.close()
            self.anal_h5file = None
        except AttributeError:
            pass
        
    def require_analysis_group(self, groupname, replace=False):
        self.open_analysis_backing()
        if replace:
            try:
                del self.anal_h5file[groupname]
            except KeyError:
                pass
        return self.anal_h5file.require_group(groupname)
    