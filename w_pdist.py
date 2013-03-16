from __future__ import print_function, division; __metaclass__ = type
import logging
from itertools import izip
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
import numpy, h5py
from fasthist import histnd, normhistnd
from westtools import h5io

log = logging.getLogger('westtools.w_pcpdist')

def isiterable(x):
    try:
        iter(x)
    except TypeError:
        return False
    else:
        return True
    
class PcoordHistHelper:
    '''Create histograms from WEST progress coordinate data. This class coordinates reading
    data from the HDF5 file and performing the histogram binning. The start and stop
    iterations must be provided, as must be a data reader instance.
    
    To use:
      1) Call ``construct_bins`` to build the bins for the histogram. This will trigger a scan
         (possibly very expensive in I/O) of all the progress coordinate data unless explicit
         boundaries are supplied. After this, bin boundaries are available in ``binbounds``
         and midpoints of bins are available in ``midpoints``.
      2) Call ``construct_histogram`` to build the histogram.
      3) The histogram as a function of iteration is then available as ``histograms``, and
         the histogram averaged over iterations is available as ``avg_histogram``.
         
    Input is taken by entire iteration when feasible, and one segment (across timepoints)
    or one timepoint (across segments) otherwise. The switch from fast (but high-memory)
    per-iteration reads to slower, smaller reads occurs when the number of entries in the
    progress coordinate data set for one iteration exceeds ``max_n_elems``, which defaults
    to 100 million (for ~400 MB of single-precision data or ~800 MB of double-precision data)
    but may be adjusted either at the class level or the instance level.
    '''
    
    max_n_elems = 100000000
        
    def __init__(self, iter_start, iter_stop, data_reader):
        self.data_reader = data_reader
        self.iter_start = iter_start
        self.iter_stop = iter_stop
        
        self.ndim = None
        self.ntimepoints = None
        self.binbounds = None  # bin boundaries for each dimension
        self.midpoints = None  # bin midpoints for each dimension 
        self.data_range = None # data range for each dimension, as the pairs (min,max)
        self.histograms = None  # Final histogram time series
        self.avg_histogram = None # Final histogram
        
    def parse_binspec(self, binspec):
        namespace = {'numpy': numpy,
                     'inf': float('inf')}
                     
        try:
            binspec_compiled = eval(binspec,namespace)
        except Exception as e:
            raise ValueError('invalid bin specification: {!r}'.format(e))
        return binspec_compiled
    
        
    def construct_bins(self, bins):
        '''
        Construct bins according to ``bins``, which may be:
        
          1) A scalar integer (for that number of bins in each dimension)
          2) A sequence of integers (specifying number of bins for each dimension)
          3) A sequence of sequences of bin boundaries (specifying boundaries for each dimension)
          
        Sets ``self.binbounds`` to a list of arrays of bin boundaries appropriate for passing to 
        fasthist.histnd, along with ``self.midpoints`` to the midpoints of the bins.
        '''
        
        if not isiterable(bins):
            self._construct_bins_from_scalar(bins)
        elif not isiterable(bins[0]):
            self._construct_bins_from_int_seq(bins)
        else:
            self._construct_bins_from_bound_seqs(bins)
            
        if log.isEnabledFor(logging.DEBUG):
            log.debug('binbounds: {!r}'.format(self.binbounds))
            
    def scan_data_shape(self):
        if self.ndim is None:
            pcoord_ds = self.data_reader.get_iter_group(self.iter_start)['pcoord']
            self.ntimepoints = pcoord_ds.shape[1]
            self.ndim = pcoord_ds.shape[2]
            self.pcoord_dtype = pcoord_ds.dtype
        
            
    def scan_data_range(self):
        '''Scan input data for range in each dimension. The number of dimensions is determined
        from the shape of the progress coordinate as of self.iter_start.'''
        
        self.scan_data_shape()
        
        pcoord_ds = self.data_reader.get_iter_group(self.iter_start)['pcoord']
        self.ntimepoints = pcoord_ds.shape[1]
        ndim = self.ndim = pcoord_ds.shape[2]
        pcoord_dtype = pcoord_ds.dtype
        
        try:
            minval = numpy.finfo(pcoord_dtype).min
            maxval = numpy.finfo(pcoord_dtype).max
        except ValueError:
            minval = numpy.iinfo(pcoord_dtype).min
            maxval = numpy.iinfo(pcoord_dtype).max
        
        data_range = self.data_range = [tuple((maxval,minval)) for _i in xrange(ndim)]
        
        log.info('determining minimum and maximum values for {} dimensions across {} iterations'
                 .format(ndim,(self.iter_stop - self.iter_start)))
        for n_iter in xrange(self.iter_start, self.iter_stop):
            if log.isEnabledFor(logging.DEBUG):
                log.debug('scanning iteration {}'.format(n_iter))
            pcoord_ds = self.data_reader.get_iter_group(n_iter)['pcoord']
            if numpy.multiply.reduce(pcoord_ds.shape) > self.max_n_elems:
                for seg_id in xrange(pcoord_ds.shape[0]):
                    for idim in xrange(ndim):
                        dimdata = pcoord_ds[seg_id,:,idim]
                        current_min, current_max = data_range[idim]
                        current_min = min(current_min, dimdata.min())
                        current_max = max(current_max, dimdata.max())
                        data_range[idim] = (current_min, current_max)
                        del dimdata
            else:
                for idim in xrange(ndim):
                    dimdata = pcoord_ds[:,:,idim]
                    current_min, current_max = data_range[idim]
                    current_min = min(current_min, dimdata.min())
                    current_max = max(current_max, dimdata.max())
                    data_range[idim] = (current_min, current_max)
                    del dimdata
            del pcoord_ds
        log.debug('data ranges: {!r}'.format(data_range))
                
    def _construct_bins_from_scalar(self, bins):
        if self.data_range is None:
            self.scan_data_range()        

        self.binbounds = []
        self.midpoints = []        
        for idim in xrange(self.ndim):
            lb, ub = self.data_range[idim]
            # Advance just beyond the upper bound of the range, so that we catch 
            # the maximum in the histogram
            ub *= 1.01
            
            boundset = numpy.linspace(lb,ub,bins+1)
            midpoints = (boundset[:-1] + boundset[1:]) / 2.0
            self.binbounds.append(boundset)
            self.midpoints.append(midpoints)
            
    def _construct_bins_from_int_seq(self, bins):
        if self.data_range is None:
            self.scan_data_range()        

        self.binbounds = []
        self.midpoints = []        
        for idim in xrange(self.ndim):
            lb, ub = self.data_range[idim]
            # Advance just beyond the upper bound of the range, so that we catch 
            # the maximum in the histogram
            ub *= 1.01
            
            boundset = numpy.linspace(lb,ub,bins[idim]+1)
            midpoints = (boundset[:-1] + boundset[1:]) / 2.0
            self.binbounds.append(boundset)
            self.midpoints.append(midpoints)
               
    def _construct_bins_from_bound_seqs(self, bins):
        self.binbounds = []
        self.midpoints = []
        for boundset in bins:
            boundset = numpy.asarray(boundset)
            if (numpy.diff(boundset) <= 0).any():
                raise ValueError('boundary set {!r} is not strictly monotonically increasing'.format(boundset))
            self.binbounds.append(boundset)
            self.midpoints.append((boundset[:-1]+boundset[1:])/2.0)
            
    def construct_histogram(self):
        '''Construct a histogram using bins previously constructed with ``construct_bins()``.
        The time series of histogram values is stored in ``histograms`` and the average over
        time is stored in ``avg_histogram``. Each histogram in the time series is normalized,
        as is the average histogram.'''
        
        self.scan_data_shape()
        
        iter_count = self.iter_stop - self.iter_start
        histograms = self.histograms = numpy.zeros((iter_count,) + tuple(len(binbounds)-1 for binbounds in self.binbounds),
                                                 dtype=numpy.float64)
        binbounds = [numpy.require(boundset, self.pcoord_dtype, 'C') for boundset in self.binbounds]
        for iiter, n_iter in enumerate(xrange(self.iter_start, self.iter_stop)):
            if log.isEnabledFor(logging.DEBUG):
                log.debug('binning iteration {}'.format(n_iter))
            iter_group = self.data_reader.get_iter_group(n_iter)
            pcoord_ds = iter_group['pcoord']
            
            npts = pcoord_ds.shape[1]
            weights = iter_group['seg_index']['weight']
            initpoint = 1 if n_iter != self.iter_start else 0
            
            # If we have too much data, then loop over it rather than load it all at once
            # The default max_n_elems is 100 million, or ~800 megs of double-precision floats
            # HDF5 reads are the slowest part of this procedure, so we minimize them by 
            # assuming that the number of timepoints is less than the number of segments, and
            # slicing by timepoint rather than segment.
            if numpy.multiply.reduce(pcoord_ds.shape) > self.max_n_elems:
                for ipt in xrange(initpoint,npts): 
                    pcoord_data = pcoord_ds[:,ipt,:]
                    histnd(pcoord_data, binbounds, weights, out=histograms[iiter], binbound_check = False)
                    del pcoord_data
            else:
                pcoord_data = pcoord_ds[:,initpoint:,:]
                for ipt in xrange(npts-initpoint):
                    histnd(pcoord_data[:,ipt,:], binbounds, weights, out=histograms[iiter], binbound_check = False)
                del pcoord_data
            
            del iter_group, weights, pcoord_ds
            
            # normalize histogram
            normhistnd(histograms[iiter],binbounds)
            
        # Calculate average histogram and normalize
        # This is fast, and may be convenient for quick visualization.
        self.avg_histogram = histograms.sum(axis=0)
        normhistnd(self.avg_histogram, binbounds)
        

class WPCPDist(WESTTool):
    prog='w_pcpdist'
    description = '''\
Calculate probability distribution of progress coordinate values and its time evolution. 
'''
    
    def __init__(self):
        super(WPCPDist,self).__init__()
        
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection(self.data_reader)
        self.iter_range.include_args['iter_step'] = False
        
        self.binspec = None
        self.output_file = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
                
        parser.add_argument('-b', '--bins', dest='bins', metavar='BINEXPR', default='100',
                            help='''Use BINEXPR for bins. This may be an integer, which will be used for each
                            dimension of the progress coordinate; a list of integers (formatted as [n1,n2,...])
                            which will use n1 bins for the first dimension, n2 for the second dimension, and so on;
                            or a list of lists of boundaries (formatted as [[a1, a2, ...], [b1, b2, ...], ... ]), which
                            will use [a1, a2, ...] as bin boundaries for the first dimension, [b1, b2, ...] as bin boundaries
                            for the second dimension, and so on. (Default: 100 bins in each dimension.)''')
        
        parser.add_argument('-o', '--output', dest='output', default='pcpdist.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')



    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open()
        self.iter_range.process_args(args)
        self.binspec = args.bins        
        self.output_file = h5py.File(args.output, 'w')
        h5io.stamp_creator_data(self.output_file)
    
    def go(self):
        hh = PcoordHistHelper(self.iter_range.iter_start, self.iter_range.iter_stop, self.data_reader)
        hh.construct_bins(hh.parse_binspec(self.binspec))
        hh.construct_histogram()
        
        for idim, (binbounds, midpoints) in enumerate(izip(hh.binbounds, hh.midpoints)):
            self.output_file['binbounds_{}'.format(idim)] = binbounds
            self.output_file['midpoints_{}'.format(idim)] = midpoints
            
        self.output_file['histograms'] = hh.histograms
        self.output_file['avg_histogram'] = hh.avg_histogram
        self.output_file['n_iter'] = self.iter_range.iter_range()
        self.iter_range.record_data_iter_range(self.output_file['avg_histogram'])
        self.iter_range.record_data_iter_range(self.output_file['histograms'])
        self.output_file.close()

if __name__ == '__main__':
    WPCPDist().main()
    
