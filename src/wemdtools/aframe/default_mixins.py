from __future__ import division, print_function; __metaclass__ = type

import logging

log = logging.getLogger(__name__)

import wemd, wemdtools

from wemdtools.aframe import AnalysisMixin, ArgumentError
                    
class IterRangeMixin(AnalysisMixin):
    '''A mixin for limiting the range of data considered for a given analysis. This should go after
    DataManagerMixin'''
    def __init__(self):
        super(IterRangeMixin,self).__init__()
        
        self.first_iter = None
        self.last_iter = None
        self.iter_step = None

    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc = super(IterRangeMixin,self).add_args
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
    
    def process_args(self, args, upcall = True):
        self.first_iter = args.first_iter or 1
        self.last_iter = args.last_iter
        self.iter_step = args.iter_step or 1

        if upcall:
            try:
                upfunc = super(IterRangeMixin,self).process_args
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
        
    def iter_block_iter(self):
        '''Return an iterable of (block_begin,past_block_end) over the blocks of iterations
        selected by --first/--last/--step.'''
                            
        for blkfirst in xrange(self.first_iter, self.last_iter+1, self.iter_step):
            yield(blkfirst, min(self.last_iter, blkfirst+self.iter_step-1)+1)
             
        
    def n_iter_blocks(self):
        npoints = self.last_iter - self.first_iter + 1
        if npoints % self.iter_step == 0:
            return npoints // self.iter_step
        else:
            return npoints // self.iter_step + 1
