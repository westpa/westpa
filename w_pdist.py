from __future__ import print_function, division; __metaclass__ = type
import sys, os, logging
from wt2.tool_classes import WEMDTool, WEMDDataReader, SegSelector
import numpy, h5py
from fasthist import hist


import wemd
from wemd.data_manager import (weight_dtype, n_iter_dtype, seg_id_dtype, utime_dtype, vstr_dtype, 
                               istate_type_dtype, istate_status_dtype)

log = logging.getLogger('wemdtools.' + __name__)

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
            log.debug('finding minimum and maximum values')
            minbound, maxbound = self.get_min_max()
            self.binbounds = numpy.linspace(minbound, maxbound, bins+1)
            log.info('using {:d} bins on [{:g},{:g}]'.format(bins, minbound, maxbound))
        else:
            self.binbounds = numpy.asarray(bins)
            log.info('using {:d} manually-specified bins'.format(len(self.binbounds)-1))
            
        self.midpoints = (self.binbounds[:-1] + self.binbounds[1:]) / 2.0
        self.nbins = len(self.binbounds)-1
        assert len(self.midpoints) ==  self.nbins
        self.hist = None
        self.normfactor = None
            
    def get_min_max(self):
        '''Scan input data for minimum and maximum values'''
        
        # Yes, this is correct; we want minval to be guaranteed to decrease with any
        # valid data point, and so on.
        minval = numpy.finfo(numpy.float64).max
        maxval = numpy.finfo(numpy.float64).min
        
        for n_iter in xrange(self.segment_selection.start_iter, self.segment_selection.stop_iter):
            
            if self.dsspec.indexed:
                for seg_id in self.segment_selection.from_iter(n_iter):
                    val = self.dsspec[n_iter,seg_id]
                    minval = min(minval, val.min())
                    maxval = max(maxval, val.max())
            else:
                seg_ids = self.segment_selection.from_iter(n_iter)
                if seg_ids:
                    values = numpy.array([self.dsspec[n_iter,seg_id] for seg_id in seg_ids])
                    minval = min(minval,values.min())
                    maxval = max(maxval,values.max())
                
                
                # Retrieve entire iteration's values
                #values = self.dsspec[n_iter]
                
                # But only consider those segments we've been asked to
                #for seg_id in self.segment_selection.from_iter(n_iter):
                #    minval = min(minval,values[seg_id].min())
                #    maxval = max(maxval,values[seg_id].max())
                    del values
                del seg_ids
        
        return minval, maxval
    
    def construct_histogram(self):
        whist = self.hist = numpy.zeros((self.nbins,), numpy.float64)
        
        dsspec = self.dsspec
        binbounds = self.binbounds
        dx = numpy.diff(binbounds)
        if (dx <= 0.0).any():
            raise ValueError('bin bounds are not strictly monotonic')
        
        start_iter = self.segment_selection.start_iter
        stop_iter  = self.segment_selection.stop_iter
        
        log.info('constructing histogram over {1:d} segments spanning {0:d} iterations'
                 .format(stop_iter-start_iter, len(self.segment_selection)))
        
        for n_iter in xrange(start_iter, stop_iter):
            seg_ids = sorted(self.segment_selection.from_iter(n_iter))
            weights = self.data_manager.get_weights(n_iter, seg_ids)
            if dsspec.indexed:
                for iseg, seg_id in enumerate(seg_ids):
                    values = numpy.asarray(dsspec[n_iter,seg_id], dtype=numpy.float64)
                    try:
                        hist(values, binbounds, out=whist, cweight=weights[iseg], binbound_check=False)
                    except:
                        log.exception('exception for segment {}:{}'.format(n_iter,seg_id))
                        raise
                    del values
            else:
                all_values = dsspec[n_iter]
                for iseg, seg_id in enumerate(seg_ids):
                    values = numpy.asarray(all_values[seg_id], dtype=numpy.float64)
                    try:
                        hist(values, binbounds, out=whist, cweight=weights[iseg], binbound_check=False)
                    except:
                        log.exception('exception for segment {}:{}'.format(n_iter,seg_id))
                        raise
                    del values
                del all_values
            
            del seg_ids, weights

        I = (whist * dx).sum()
        self.normfactor = 1.0/I
        whist /= I
        

class WPDistTool(WEMDTool):
    prog='w_pdist'
    description = '''\
Calculate probability distributions over a weighted ensemble data set. 
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
        hh = HistogramHelper(self.dsspec, self.bins, self.segselector.segment_selection, 
                             self.data_reader.data_manager)
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
    
