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
from oldtools.aframe import AnalysisMixin

class BinningMixin(AnalysisMixin):
    '''A mixin for performing binning on WEST data.'''
    
    def __init__(self):
        super(BinningMixin,self).__init__()
        
        self.mapper = None
        self.n_bins = None
        
        self.discard_bin_assignments = False
        self.binning_h5gname = 'binning'
        self.binning_h5group = None
        self.mapper_hash = None

    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc = super(BinningMixin,self).add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
        
        group = parser.add_argument_group('binning options')
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--binexpr', '--binbounds', dest='binexpr',
                            help='''Construct rectilinear bins from BINEXPR. This must be a list of lists of bin boundaries
                            (one list of bin boundaries for each dimension of the progress coordinate), formatted as a Python 
                            expression. E.g. "[[0,1,2,4,inf], [-inf,0,inf]]".''')
        group.add_argument('--discard-bin-assignments', dest='discard_bin_assignments', action='store_true',
                           help='''Discard any existing bin assignments stored in the analysis HDF5 file.''')
    
    def process_args(self, args, upcall = True):        
        if args.binexpr:
            westpa.rc.pstatus("Constructing rectilinear bin boundaries from the following expression: '{}'".format(args.binexpr))
            self.mapper = self.mapper_from_expr(args.binexpr)
        else:
            westpa.rc.pstatus('Loading bin boundaries from WEST system')
            system = westpa.rc.get_system_driver()
            self.mapper = system.bin_mapper
            
        self.n_bins = self.mapper.nbins
        _pdat, self.mapper_hash = self.mapper.pickle_and_hash()
        westpa.rc.pstatus('  {:d} bins'.format(self.n_bins))
        westpa.rc.pstatus('  identity hash {}'.format(self.mapper_hash))
        
        self.discard_bin_assignments = bool(args.discard_bin_assignments)
        
        if upcall:
            try:
                upfunc = super(BinningMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)
        
    def mapper_from_expr(self, expr):
        from westpa.binning import RectilinearBinMapper
        namespace = {'numpy': numpy,
                     'inf': float('inf')}
        
        try:
            return RectilinearBinMapper(eval(expr, namespace))
        except TypeError as e:
            if 'has no len' in str(e):
                raise ValueError('invalid bin boundary specification; you probably forgot to make a list of lists')

    def write_bin_labels(self, dest, 
                         header='# bin labels:\n', 
                         format='# bin {bin_index:{max_iwidth}d} -- {label!s}\n'):
        '''Print labels for all bins in ``self.mapper`` to ``dest``.  If provided, ``header`` 
        is printed before any labels.   The ``format`` string specifies how bin labels are to be printed.  Valid entries are:
          * ``bin_index`` -- the zero-based index of the bin
          * ``label`` -- the label, as obtained by ``bin.label``
          * ``max_iwidth`` -- the maximum width (in characters) of the bin index, for pretty alignment
        '''
        dest.write(header or '')
        max_iwidth = len(str(self.mapper.nbins-1))
        for (ibin, label) in enumerate(self.mapper.labels):
            dest.write(format.format(bin_index=ibin, label=label, max_iwidth=max_iwidth))
    
    def require_binning_group(self):
        if self.binning_h5group is None:
            self.binning_h5group = self.anal_h5file.require_group(self.binning_h5gname)
        return self.binning_h5group
    
    def delete_binning_group(self):
        self.binning_h5group = None
        del self.anal_h5file[self.binning_h5gname]

    def record_data_binhash(self, h5object):
        '''Record the identity hash for self.mapper as an attribute on the given HDF5 object (group or dataset)'''
        h5object.attrs['binhash'] = self.mapper_hash
        
    def check_data_binhash(self, h5object):
        '''Check whether the recorded bin identity hash on the given HDF5 object matches the identity hash for self.mapper'''
        return h5object.attrs.get('binhash') == self.mapper_hash 
            
    def assign_to_bins(self):
        '''Assign WEST segment data to bins.  Requires the DataReader mixin to be in the inheritance tree'''
        self.require_binning_group()        
        
        n_iters = self.last_iter - self.first_iter + 1
        max_n_segs = self.max_iter_segs_in_range(self.first_iter, self.last_iter)
        pcoord_len = self.get_pcoord_len(self.first_iter)
        
        assignments = numpy.zeros((n_iters, max_n_segs,pcoord_len), numpy.min_scalar_type(self.n_bins))
        populations = numpy.zeros((n_iters, pcoord_len, self.n_bins), numpy.float64)
        
        westpa.rc.pstatus('Assigning to bins...')
        
        for (iiter, n_iter) in enumerate(range(self.first_iter, self.last_iter+1)):
            westpa.rc.pstatus('\r  Iteration {:d}'.format(n_iter), end='')
            seg_index = self.get_seg_index(n_iter)
            pcoords = self.get_iter_group(n_iter)['pcoord'][...]
            weights = seg_index['weight']
            
            for seg_id in range(len(seg_index)):
                assignments[iiter,seg_id,:] = self.mapper.assign(pcoords[seg_id,:,:])
            
            for it in range(pcoord_len):
                populations[iiter, it, :] = numpy.bincount(assignments[iiter,:len(seg_index),it], weights, minlength=self.n_bins)
        
            westpa.rc.pflush()
            del pcoords, weights, seg_index
         
        assignments_ds = self.binning_h5group.create_dataset('bin_assignments', data=assignments, compression='gzip')
        populations_ds = self.binning_h5group.create_dataset('bin_populations', data=populations, compression='gzip')
        
        for h5object in (self.binning_h5group, assignments_ds, populations_ds):
            self.record_data_iter_range(h5object)
            self.record_data_iter_step(h5object, 1)
            self.record_data_binhash(h5object)
                
        westpa.rc.pstatus()
            
    def require_bin_assignments(self):
        self.require_binning_group()
        do_assign = False
        if self.discard_bin_assignments:
            westpa.rc.pstatus('Discarding existing bin assignments.')
            do_assign = True
        elif 'bin_assignments' not in self.binning_h5group:
            do_assign = True
        elif not self.check_data_iter_range_least(self.binning_h5group):
            westpa.rc.pstatus('Existing bin assignments are for incompatible first/last iterations; deleting assignments.')
            do_assign = True
        elif not self.check_data_binhash(self.binning_h5group):
            westpa.rc.pstatus('Bin definitions have changed; deleting existing bin assignments.')
            do_assign = True
    
        if do_assign:
            self.delete_binning_group()
            self.assign_to_bins()
        else:
            westpa.rc.pstatus('Using existing bin assignments.')
            
    def get_bin_assignments(self, first_iter = None, last_iter = None):
        return self.slice_per_iter_data(self.binning_h5group['bin_assignments'], first_iter, last_iter)

    def get_bin_populations(self, first_iter = None, last_iter = None):
        return self.slice_per_iter_data(self.binning_h5group['bin_populations'], first_iter, last_iter)
    