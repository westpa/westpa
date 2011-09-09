from __future__ import division, print_function; __metaclass__ = type

import logging

log = logging.getLogger(__name__)

import numpy

import wemd
from wemdtools.aframe import AnalysisMixin

class BinningMixin(AnalysisMixin):
    '''A mixin for performing binning on WEMD data.
    
    Requires: IterRangeMixin'''
    
    def __init__(self):
        super(BinningMixin,self).__init__()
        
        self.region_set = None
        self.n_bins = None
        self.n_dim = None
        
        self.__discard_bin_assignments = False
        self.__bins_checked = False
        self.binning_h5gname = 'binning'
        self.binning_h5group = None
        self.region_set_hash = None

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
            wemd.rc.pstatus("Constructing rectilinear bin boundaries from the following expression: '{}'".format(args.binexpr))
            self.region_set = self.region_set_from_expr(args.binexpr)
        else:
            wemd.rc.pstatus('Loading bin boundaries from WEMD system')
            system = wemd.rc.get_system_driver()
            self.region_set = system.new_region_set()
            
        self.n_bins = len(self.region_set.get_all_bins())
        self.n_dim = self.region_set.n_dim
        self.region_set_hash = self.region_set.identity_hash()
        wemd.rc.pstatus('  {:d} bins in {:d} dimension(s)'.format(self.n_bins, self.n_dim))
        wemd.rc.pstatus('  identity hash {}'.format(self.region_set_hash.hexdigest()))
        
        self.__discard_bin_assignments = bool(args.discard_bin_assignments)
        
        if upcall:
            try:
                upfunc = super(BinningMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)
        
    def region_set_from_expr(self, expr):
        from wemd.pcoords import RectilinearRegionSet, PiecewiseRegionSet
        namespace = {'numpy': numpy,
                     'RectilinearRegionSet': RectilinearRegionSet,
                     'PiecewiseRegionSet': PiecewiseRegionSet,
                     'inf': float('inf')}
        
        return RectilinearRegionSet(eval(expr, namespace))
        
    def write_bin_labels(self, dest, 
                         header='# bin labels:\n', 
                         format='# bin {bin_index:{max_iwidth}d} -- {label!s}\n'):
        '''Print labels for all bins in the given RegionSet (or ``self.region_set``) to ``dest``.  If provided, ``header`` 
        is printed before any labels.   The ``format`` string specifies how bin labels are to be printed.  Valid entries are:
          * ``bin_index`` -- the zero-based index of the bin
          * ``label`` -- the label, as obtained by ``bin.label``
          * ``max_iwidth`` -- the maximum width (in characters) of the bin index, for pretty alignment
        '''
        dest.write(header or '')
        bins = self.region_set.get_all_bins()
        max_iwidth = len(str(len(bins)-1))
        for (ibin, bin) in enumerate(bins):
            dest.write(format.format(bin_index=ibin, label=bin.label, max_iwidth=max_iwidth))
    
    def __require_group(self):
        if self.binning_h5group is None:
            self.binning_h5group = self.anal_h5file.require_group(self.binning_h5gname)
        return self.binning_h5group
    
    def __delete_group(self):
        self.binning_h5group = None
        del self.anal_h5file[self.binning_h5gname]

    def record_data_binhash(self, h5object):
        '''Record the identity hash for self.region_set as an attribute on the given HDF5 object (group or dataset)'''
        h5object.attrs['binhash'] = self.region_set_hash.digest()
        
    def check_data_binhash(self, h5object):
        '''Check whether the recorded bin identity hash on the given HDF5 object matches the identity hash for self.region_set'''
        return h5object.attrs.get('binhash') == self.region_set_hash.digest() 
                
    def check_bin_data(self):
        '''Check to see that existing binning data corresponds to the same bin topology and iteration range as requested'''
        
        self.__require_group()
        
        if self.__discard_bin_assignments:
            wemd.rc.pstatus('Discarding existing binning data.')
            self.__delete_group()
        elif 'bin_assignments' in self.binning_h5group:
            if not self.check_data_iter_range_least(self.binning_h5group):
                wemd.rc.pstatus('Existing binning data is for incompatible first/last iterations; deleting.')
                self.__delete_group()
            elif not self.check_data_binhash(self.binning_h5group):
                wemd.rc.pstatus('Bin definitions have changed; deleting existing binning data.')
                self.__delete_group()
        
        self.__require_group()
        
    def assign_to_bins(self):
        '''Requires the DataReader mixin to be in the inheritance tree'''
        self.__require_group()        
        
        n_iters = self.last_iter - self.first_iter + 1
        max_n_segs = self.max_iter_segs_in_range(self.first_iter, self.last_iter)
        pcoord_len = self.get_pcoord_len(self.first_iter)
        
        assignments = numpy.zeros((n_iters, max_n_segs,pcoord_len), numpy.min_scalar_type(self.n_bins))
        populations = numpy.zeros((n_iters, pcoord_len, self.n_bins), numpy.float64)
        
        wemd.rc.pstatus('Assigning to bins...')
        
        for (iiter, n_iter) in enumerate(xrange(self.first_iter, self.last_iter+1)):
            wemd.rc.pstatus('\r  Iteration {:d}'.format(n_iter), end='')
            seg_index = self.get_seg_index(n_iter)
            pcoords = self.get_iter_group(n_iter)['pcoord'][...]
            weights = seg_index['weight']
            
            for seg_id in xrange(len(seg_index)):
                assignments[iiter,seg_id,:] = self.region_set.map_to_all_indices(pcoords[seg_id,:,:])
            
            for it in xrange(pcoord_len):
                populations[iiter, it, :] = numpy.bincount(assignments[iiter,:len(seg_index),it], weights, minlength=self.n_bins)
        
            wemd.rc.pflush()
            del pcoords, weights, seg_index
         
        assignments_ds = self.binning_h5group.create_dataset('bin_assignments', data=assignments, compression='gzip')
        populations_ds = self.binning_h5group.create_dataset('bin_populations', data=populations, compression='gzip')
        
        for h5object in (self.binning_h5group, assignments_ds, populations_ds):
            self.record_data_iter_range(h5object)
            self.record_data_iter_step(h5object, 1)
            self.record_data_binhash(h5object)
                
        wemd.rc.pstatus()
            
    def require_bin_assignments(self):
        self.check_bin_data()
                
        # The group will be empty if the user requested the data to go away, or if the data is not conformant
        # with the current bins and iteration range, so the following lets us know if we need
        # to recalculate 
        if not ('bin_assignments' in self.binning_h5group and 'bin_populations' in self.binning_h5group):
            self.assign_to_bins()
        else:
            wemd.rc.pstatus('Using existing binning data.')
    
    def get_bin_assignments(self, first_iter = None, last_iter = None):
        return self.slice_per_iter_data(self.binning_h5group['bin_assignments'], first_iter, last_iter)

    def get_bin_populations(self, first_iter = None, last_iter = None):
        return self.slice_per_iter_data(self.binning_h5group['bin_populations'], first_iter, last_iter)
        
