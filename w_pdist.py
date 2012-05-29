from __future__ import print_function, division; __metaclass__ = type
import sys, os, logging
from wt2.tool_classes import WEMDTool, WEMDDataReader, SegSelector
import numpy, h5py
from fasthist import hist


import wemd
from wemd.data_manager import (weight_dtype, n_iter_dtype, seg_id_dtype, utime_dtype, vstr_dtype, 
                               istate_type_dtype, istate_status_dtype)

log = logging.getLogger(__name__)

def parse_binspec(binspec):
    namespace = {'inf': float('inf'),
                 'numpy': numpy,
                 'os': os,
                 'sys': sys}
    try:
        return eval(binspec, namespace)
    except Exception as e:
        raise ValueError('invalid bin boundary expression {!r}: {!r}'.format(binspec,e))
    
class HistogramHelper:
    def __init__(self, dsspec, bins, segment_selection=None, data_manager = None):
        self.data_manager = data_manager or wemd.rc.get_data_manager()
        self.dsspec = dsspec
        self.segment_selection = segment_selection
        
        if numpy.isscalar(bins):
            minbound, maxbound = self.get_min_max()
            self.binbounds = numpy.linspace(minbound, maxbound, bins)
        else:
            self.binbounds = numpy.asarray(bins)
            
        self.midpoints = (self.binbounds[:-1] + self.binbounds[1:]) / 2.0
        self.nbins = len(self.midpoints)
        self.hist = None
        self.normfactor = None
            
    def get_min_max(self):
        '''Scan input data for minimum and maximum values'''
        
        # Yes, this is correct; we want minval to be guaranteed to decrease with any
        # valid data point, and so on.
        minval = numpy.finfo(numpy.float64).max
        maxval = numpy.finfo(numpy.float64).min
        
        for n_iter in xrange(self.segment_selection.start_iter, self.segment_selection.stop_iter):
            values = self.dsspec[n_iter]
            minval = min(minval,values.min())
            maxval = max(maxval,values.max())
            del values
        
        return minval, maxval
    
    def construct_histogram(self):
        whist = self.hist = numpy.zeros((self.nbins,), numpy.float64)
        
        dsspec = self.dsspec
        binbounds = self.binbounds
        dx = numpy.diff(binbounds)
        if (dx <= 0.0).any():
            raise ValueError('bin bounds are not strictly monotonic')
        
        for n_iter in xrange(self.segment_selection.start_iter, self.segment_selection.stop_iter):
            seg_ids = sorted(self.segment_selection.from_iter(n_iter))
            weights = self.data_manager.get_weights(n_iter, seg_ids)
            for iseg, seg_id in enumerate(seg_ids):
                values = numpy.asarray(dsspec[n_iter,seg_id], dtype=numpy.float64)
                try:
                    hist(values, binbounds, out=whist, cweight=weights[iseg], binbound_check=False)
                except:
                    log.exception('exception for segment {}:{}'.format(n_iter,seg_id))
                    raise
                del values
            
            del seg_ids, weights
        
        I = (whist * dx).sum()
        self.normfactor = 1.0/I
        whist /= I
        

class WPDistTool(WEMDTool):
    prog='w_pdist'
    description = '''\
Calculate probability distributions over a weighted ensemble data set. 

Individual datasets can be selected for writing using the -d/--dataset option
(which may be specified more than once). The simplest form is ``-d dsname``,
which causes data from dataset ``dsname`` along the trace to be stored to
HDF5.  The dataset is assumed to be stored on a per-iteration basis, with
the first dimension corresponding to seg_id and the second dimension
corresponding to time within the segment.  Further options are specified
as comma-separated key=value pairs after the data set name, as in

    -d dsname;alias=newname;index=idsname;file=otherfile.h5;slice=[100,...]
    
The following options for datasets are supported:

    alias=newname
        When writing this data to HDF5 or text files, use ``newname``
        instead of ``dsname`` to identify the dataset. This is mostly of
        use in conjunction with the ``slice`` option in order, e.g., to
        retrieve two different slices of a dataset and store then with
        different names for future use.

    index=idsname
        The dataset is not stored on a per-iteration basis for all
        segments, but instead is stored as a single dataset whose
        first dimension indexes n_iter/seg_id pairs. The index to
        these n_iter/seg_id pairs is ``idsname``.
    
    file=otherfile.h5
        Instead of reading data from the main WEMD HDF5 file (usually
        ``wemd.h5``), read data from ``otherfile.h5``.
        
    slice=[100,...]
        Retrieve only the given slice from the dataset. This can be
        used to pick a subset of interest to minimize I/O.
        
-------------------------------------------------------------------------------
'''
    
    def __init__(self):
        super(WPDistTool,self).__init__()
        
        self.data_reader = WEMDDataReader()
        self.segselector = SegSelector()
        
        self.dsspec = None
        self.bins = None
        self.output_file = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.segselector.add_args(parser)
        
        parser.add_argument('dataset',
                            #this breaks argparse (see http://bugs.python.org/issue11874) 
                            #metavar='DSNAME[,alias=ALIAS][,index=INDEX][,file=FILE][,slice=SLICE]',
                            metavar='DSNAME',
                            help='''Calculate distributions for the dataset named DSNAME. An extended form like
                            DSNAME[;alias=ALIAS][;index=INDEX][;file=FILE][;slice=SLICE] will
                            obtain the dataset from the given FILE instead of the main WEMD HDF5 file,
                            slice it by SLICE, call it ALIAS in output, and/or access per-segment data by a n_iter,seg_id
                            INDEX instead of a seg_id indexed dataset in the group for iteration n_iter.''')
        
        bgroup = parser.add_mutually_exclusive_group()
        bgroup.add_argument('-N', '--nbins', dest='nhistbins', type=int, default=1000,
                            help='''Use NHISTBINS bins to bin data for histograms (default: %(default)s).''')
        bgroup.add_argument('--binexpr', dest='binexpr',
                            help='''Evaluate BINEXPR and use as bin boundaries.''')
        
        parser.add_argument('-o', '--output', dest='output', default='pdist.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open()        
        self.segselector.process_args(args)
        self.output_file = h5py.File(args.output, 'w')
        self.dsspec = self.data_reader.parse_dssel_str(args.dataset)
        
        if not args.binexpr:
            self.bins = args.nhistbins
        else:
            self.bins = parse_binspec(args.binexpr)
    
    def go(self):
        hh = HistogramHelper(self.dsspec, self.bins, self.segselector.segment_selection, self.data_reader.data_manager)
        hh.construct_histogram()        
        self.output_file['binbounds'] = hh.binbounds
        self.output_file['midpoints'] = hh.midpoints
        self.output_file['hist'] = hh.hist
        self.output_file['hist'].attrs['normfactor'] = hh.normfactor
        
        for ds in ('binbounds', 'midpoints', 'hist'):
            self.output_file[ds].attrs['source_name'] = hh.dsspec.name
            
        self.output_file.close()
        
            
        

if __name__ == '__main__':
    WPDistTool().main()
    
