from __future__ import division, print_function; __metaclass__ = type

import logging

log = logging.getLogger(__name__)

import wemd, wemdtools

from wemdtools.aframe import AnalysisMixin, ArgumentError

class DataManagerMixin(AnalysisMixin):
    '''A mixin for analysis requiring access to the HDF5 files generated during a WEMD run.'''

    def __init__(self):
        self.data_manager = None
        self.run_h5name = None
        
    def add_common_args(self, parser, upcall = True):
        if upcall:
            try:
                upcall = super(DataManagerMixin,self).add_common_args
            except AttributeError:
                pass
            else:
                upcall(parser)
        
        group = parser.add_argument_group('WEMD HDF5 options')
        group.add_argument('--no-cache', dest='no_cache_rundata', action='store_true',
                            help='''Disable all caching of WEMD run data.''')
        group.add_argument('run_h5name', nargs='?', metavar='WEMD_H5FILE',
                           help='''Take data from WEMD_H5FILE (default: read from the HDF5 file specified in wemd.cfg).''')

    def process_common_args(self, args, upcall = True):            
        if args.run_h5name:
            self.run_h5name = args.run_h5name
        else:
            wemd.rc.config.require('data.h5file')
            self.run_h5name = wemd.rc.config.get_path('data.h5file') 
        
        wemd.rc.pstatus("Using run data from '{}'".format(self.run_h5name))
        
        self.data_manager = wemd.rc.get_data_manager()
        self.data_manager.backing_file = self.run_h5name
        
        if not args.no_cache_rundata:
            log.debug('using caching data reader')
            from wemdtools.data_manager import CachingWEMDDataReader
            self.data_manager = CachingWEMDDataReader(self.data_manager)
        self.data_manager.open_backing(mode='r')
        
        if upcall:
            try:
                upfunc = super(DataManagerMixin,self).process_common_args
            except AttributeError:
                pass
            else:
                upfunc(args)

                    
class IterRangeMixin(AnalysisMixin):
    '''A mixin for limiting the range of data considered for a given analysis. This should go after
    DataManagerMixin'''
    def __init__(self):
        self.first_iter = None
        self.last_iter = None
        self.iter_step = None

    def add_common_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc = super(IterRangeMixin,self).add_common_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
        
        group = parser.add_argument_group('analysis range')
        group.add_argument('--start', '--begin', '--first', dest='first_iter', type=int, metavar='N_ITER', default=1,
                           help='''Begin analysis at iteration N_ITER (default: %(default)d).''')
        group.add_argument('--stop', '--end', '--last', dest='last_iter', type=int, metavar='N_ITER',
                           help='''Conclude analysis with N_ITER, inclusive (default: last completed iteration).''')
        group.add_argument('--step', dest='iter_step', type=int, metavar='STEP',
                           help='''Analyze/report in blocks of STEP iterations.''')
    
    def process_common_args(self, args, upcall = True):
        self.first_iter = args.first_iter or 1
        self.last_iter = args.last_iter
        self.iter_step = args.iter_step or 1

        if upcall:
            try:
                upfunc = super(IterRangeMixin,self).process_common_args
            except AttributeError:
                pass
            else:
                upfunc(args)
            
    def check_iter_range(self):
        assert hasattr(self, 'data_manager') and self.data_manager is not None
        
        self.first_iter = max(self.first_iter, 1)
        if self.last_iter is None or self.last_iter > self.data_manager.current_iteration - 1:
            self.last_iter = self.data_manager.current_iteration - 1
            
        if self.first_iter == self.last_iter:
            raise ArgumentError('first and last iterations are the same')

        wemd.rc.pstatus('Processing iterations from {self.first_iter:d} to {self.last_iter:d}, inclusive (step size {self.iter_step:d})'.format(self=self))
        
            
             
