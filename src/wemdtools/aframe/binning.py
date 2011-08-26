from __future__ import division, print_function; __metaclass__ = type

import logging

log = logging.getLogger(__name__)

import wemd
from wemdtools.aframe import AnalysisMixin

class BinnerMixin(AnalysisMixin):
    '''A mixin for performing binning on WEMD data.'''
    def __init__(self):
        self.region_set = None
        self.n_bins = None
        self.n_dim = None
        self.__discard_bin_assignments = False
        self.__h5gname = 'binning'
        self.__h5group = None

    def add_common_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc = super(BinnerMixin,self).add_common_args
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
    
    def process_common_args(self, args, upcall = True):        
        self.__discard_bin_assignments = bool(args.discard_bin_assignments)
        
        if args.binexpr:
            wemd.rc.pstatus("Constructing rectilinear bin boundaries from the following expression: '{}'".format(args.binexpr))
            self.region_set = self.region_set_from_expr(args.binexpr)
        else:
            wemd.rc.pstatus('Loading bin boundaries from WEMD system')
            system = wemd.rc.get_system_driver()
            self.region_set = system.new_region_set()
            
        self.n_bins = len(self.region_set.get_all_bins())
        self.n_dim = self.region_set.n_dim
        wemd.rc.pstatus('  {:d} bins in {:d} dimension(s)'.format(self.n_bins, self.n_dim))
        wemd.rc.pstatus('  identity hash {}'.format(self.region_set.identity_hash().hexdigest()))

        if upcall:
            try:
                upfunc = super(BinnerMixin,self).process_common_args
            except AttributeError:
                pass
            else:
                upfunc(args)
        
    def region_set_from_expr(self, expr):
        import numpy
        from wemd.pcoords import RectilinearRegionSet, PiecewiseRegionSet
        namespace = {'numpy': numpy,
                     'RectilinearRegionSet': RectilinearRegionSet,
                     'PiecewiseRegionSet': PiecewiseRegionSet,
                     'inf': float('inf')}
        
        return RectilinearRegionSet(eval(expr, namespace))
        
    def write_bin_labels(self, dest, region_set = None, 
                         header='# bin labels:\n', 
                         format='# bin {bin_index:{max_iwidth}d -- {label!s}\n'):
        '''Print labels for all bins in the given RegionSet (or ``self.region_set``) to ``dest``.  If provided, ``header`` 
        is printed before any labels.   The ``format`` string specifies how bin labels are to be printed.  Valid entries are:
          * ``bin_index`` -- the zero-based index of the bin
          * ``label`` -- the label, as obtained by ``bin.label``
          * ``max_iwidth`` -- the maximum width (in characters) of the bin index, for pretty alignment
        '''
        region_set = region_set or self.region_set
        dest.write(header or '')
        bins = region_set.get_all_bins()
        max_iwidth = len(str(len(bins)-1))
        for (ibin, bin) in enumerate(bins):
            dest.write(format.format(bin_index=ibin, label=bin.label, max_iwidth=max_iwidth))
    
    def __require_group(self):
        if self.__h5group is None:
            if self.__discard_bin_assignments:
                try:
                    del self.h5file[self.__h5gname]
                except AttributeError:
                    pass
            self.__h5group = self.h5file.require_group(self.__h5gname)
        return self.__h5group
                    
    
        