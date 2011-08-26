'''
Created on Aug 25, 2011

@author: mzwier
'''
from __future__ import division, print_function; __metaclass__ = type

import logging

log = logging.getLogger(__name__)

import os, sys
import h5py
import wemd, wemdtools

from wemdtools.aframe import AnalysisMixin, ArgumentError

class WEMDAnalysisTool:
    def __init__(self):
        super(WEMDAnalysisTool,self).__init__()
        self.config_required = False
        self.anal_h5name = None
        self.anal_h5file = None

    def add_common_args(self, parser, upcall = True):
        '''Add arguments to a parser common to all analyses of this type.'''
        if upcall:
            try:
                upfunc = super(WEMDAnalysisTool,self).add_common_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
                
        group = parser.add_argument_group('general analysis options')
        group.add_argument('-A', '--analysis-file', dest='anal_h5name', metavar='H5FILE', default='analysis.h5',
                            help='Store intermediate and final results in H5FILE (default: %(default)s).')
        
    def process_common_args(self, args, upcall = True):
        self.anal_h5name = args.anal_h5name
        
        if upcall:
            try:
                upfunc = super(WEMDAnalysisTool,self).process_common_args
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
