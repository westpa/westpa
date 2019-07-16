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

import numpy
import westpa

from oldtools.aframe import AnalysisMixin, ArgumentError

class IterRangeMixin(AnalysisMixin):
    '''A mixin for limiting the range of data considered for a given analysis. This should go after
    DataManagerMixin'''
    def __init__(self):
        super(IterRangeMixin,self).__init__()
        
        self.first_iter = None
        self.last_iter = None
        self.iter_step = 1

        include_args = self.include_args.setdefault('IterRangeMixin',{})
        include_args.setdefault('first_iter', True)
        include_args.setdefault('last_iter', True)
        include_args.setdefault('iter_step',True)


    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc = super(IterRangeMixin,self).add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
        
        group = parser.add_argument_group('analysis range')
        if self.include_args['IterRangeMixin']['first_iter']:
            group.add_argument('--start', '--begin', '--first', dest='first_iter', type=int, metavar='N_ITER', default=1,
                               help='''Begin analysis at iteration N_ITER (default: %(default)d).''')
        if self.include_args['IterRangeMixin']['last_iter']:
            group.add_argument('--stop', '--end', '--last', dest='last_iter', type=int, metavar='N_ITER',
                               help='''Conclude analysis with N_ITER, inclusive (default: last completed iteration).''')
        if self.include_args['IterRangeMixin']['iter_step']:            
            group.add_argument('--step', dest='iter_step', type=int, metavar='STEP',
                               help='''Analyze/report in blocks of STEP iterations.''')
    
    def process_args(self, args, upcall = True):
        if self.include_args['IterRangeMixin']['first_iter']:
            self.first_iter = args.first_iter or 1
        if self.include_args['IterRangeMixin']['last_iter']:
            self.last_iter = args.last_iter
        if self.include_args['IterRangeMixin']['iter_step']:
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
        
        self.first_iter = int(max(self.first_iter, 1))
        if self.last_iter is None or self.last_iter > self.data_manager.current_iteration - 1:
            self.last_iter = int(self.data_manager.current_iteration - 1)
            
        if self.first_iter == self.last_iter:
            raise ArgumentError('first and last iterations are the same')

        westpa.rc.pstatus('Processing iterations from {self.first_iter:d} to {self.last_iter:d}, inclusive (step size {self.iter_step:d})'.format(self=self))
        
    def iter_block_iter(self):
        '''Return an iterable of (block_first,block_last+1) over the blocks of iterations
        selected by --first/--last/--step.  NOTE WELL that the second of the pair follows Python
        iterator conventions and returns one past the last element of the block.'''
                            
        for blkfirst in range(self.first_iter, self.last_iter+1, self.iter_step):
            yield(blkfirst, min(self.last_iter, blkfirst+self.iter_step-1)+1)
             
        
    def n_iter_blocks(self):
        '''Return the number of blocks of iterations (as returned by ``iter_block_iter``) selected by --first/--last/--step.'''
        npoints = self.last_iter - self.first_iter + 1
        if npoints % self.iter_step == 0:
            return npoints // self.iter_step
        else:
            return npoints // self.iter_step + 1
            
    def record_data_iter_range(self, h5object, first_iter = None, last_iter = None):
        '''Store attributes ``first_iter`` and ``last_iter`` on the given HDF5 object (group/dataset)'''
        first_iter = first_iter or self.first_iter
        last_iter = last_iter or self.last_iter
        h5object.attrs['first_iter'] = first_iter
        h5object.attrs['last_iter'] = last_iter
        
    def record_data_iter_step(self, h5object, iter_step = None):
        '''Store attribute ``iter_step`` on the given HDF5 object (group/dataset).'''
        iter_step = iter_step or self.iter_step
        h5object.attrs['iter_step'] = iter_step
        
    def check_data_iter_range_least(self, h5object, first_iter = None, last_iter = None):
        '''Check that the given HDF5 object contains (as denoted by its ``first_iter``/``last_iter`` attributes) at least the
        data range specified.'''
        first_iter = first_iter or self.first_iter
        last_iter = last_iter or self.last_iter
        
        obj_first_iter = h5object.attrs.get('first_iter')
        obj_last_iter  = h5object.attrs.get('last_iter')
        
        return (obj_first_iter <= first_iter and obj_last_iter >= last_iter)
        
    def check_data_iter_range_equal(self, h5object, first_iter = None, last_iter = None):
        '''Check that the given HDF5 object contains per-iteration data for exactly the specified iterations (as denoted by the
        object's ``first_iter`` and ``last_iter`` attributes'''

        first_iter = first_iter or self.first_iter
        last_iter = last_iter or self.last_iter
        
        obj_first_iter = h5object.attrs.get('first_iter')
        obj_last_iter  = h5object.attrs.get('last_iter')
        
        return (obj_first_iter == first_iter and obj_last_iter == last_iter)
    
    def check_data_iter_step_conformant(self, h5object, iter_step = None):
        '''Check that the given HDF5 object contains per-iteration data at an iteration stride suitable for extracting data
        with the given stride.  (In other words, is the given ``iter_step`` a multiple of the stride with 
        which data was recorded.)'''
        
        iter_step = iter_step or self.iter_step
        obj_iter_step = h5object.attrs.get('iter_step')
        return (obj_iter_step % iter_step == 0)
    
    def check_data_iter_step_equal(self, h5object, iter_step = None):
        '''Check that the given HDF5 object contains per-iteration data at an iteration stride the same as
        that specified.'''
        iter_step = iter_step or self.iter_step
        obj_iter_step = h5object.attrs.get('iter_step')
        return (obj_iter_step == iter_step)
        
    def slice_per_iter_data(self, dataset, first_iter = None, last_iter = None, iter_step = None, axis=0):
        '''Return the subset of the given dataset corresponding to the given iteration range and stride. Unless
        otherwise specified, the first dimension of the dataset is the one sliced.'''
        
        first_iter = first_iter or self.first_iter
        last_iter = last_iter or self.last_iter
        iter_step = iter_step or self.iter_step
         
        ds_first_iter = dataset.attrs['first_iter']
        ds_last_iter  = dataset.attrs['last_iter']
        ds_iter_step  = dataset.attrs.get('iter_step', 1)
        
        if first_iter < ds_first_iter or last_iter > ds_last_iter or ds_iter_step % iter_step > 0:
            raise IndexError(('Cannot slice requested iterations [{:d},{:d}] (stride={:d}) from dataset {!r}'
                              +'with range [{:d},{:d}] (stride={:d}).'.format(first_iter,last_iter,iter_step,
                                                                              ds_first_iter,ds_last_iter,ds_iter_step)))
        
        dimslices = []
        for idim in range(len(dataset.shape)):
            if idim == axis:
                dimslices.append(slice(first_iter - ds_first_iter, last_iter - ds_first_iter + iter_step, iter_step))
            else:
                dimslices.append(slice(None,None,None))
        
        dimslices = tuple(dimslices)
        log.debug('slicing {!r} with {!r}'.format(dataset, dimslices))
        data = dataset[dimslices]
        log.debug('resulting data is of shape {!r}'.format(data.shape))
        return data        
        
    
    def iter_range(self, first_iter = None, last_iter = None, iter_step = None):
        first_iter = first_iter or self.first_iter
        last_iter = last_iter or self.last_iter
        iter_step = iter_step or self.iter_step
        
        return numpy.arange(first_iter, last_iter + 1, iter_step)
        
        
        