from __future__ import print_function, division; __metaclass__ = type
import logging
from itertools import izip
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
import numpy, h5py
from fasthist import histnd, normhistnd
from westtools import h5io
from westpa.extloader import get_object

log = logging.getLogger('westtools.w_pdist')

def isiterable(x):
    try:
        iter(x)
    except TypeError:
        return False
    else:
        return True
    
class DSSpec:
    def get_iter_data(self, n_iter):
        raise NotImplementedError
    
    def get_segment_data(self, n_iter, seg_id):
        return self.get_iter_data(n_iter)[seg_id]
    
class SingleDSSpec(DSSpec):
    @classmethod
    def from_string(cls, dsspec_string, default_h5file):
        dsname = None
        slice = None
        alias = None
        
        fields = dsspec_string.split(',')
        dsname = fields[0]
        
        for field in (field.strip() for field in fields[1:]):
            k,v = field.split('=')
            k = k.lower()
            if k == 'alias':
                alias = v
            elif k == 'slice':
                try:
                    slice = eval('numpy.index_exp' + v)
                except SyntaxError:
                    raise SyntaxError('invalid index expression {!r}'.format(v))
            else:
                raise ValueError('invalid dataset option {!r}'.format(k))
            
        return cls(default_h5file, dsname, alias, slice)
        
        
        
    
    def __init__(self, h5file, dsname, alias=None, slice=None):
        self.h5file = h5file
        self.dsname = dsname
        self.alias = alias or dsname
        self.slice = numpy.index_exp[slice] if slice else None
        
    def get_iter_data(self, n_iter):
        if self.slice:
            return self.h5file.get_iter_group(n_iter)[self.dsname][numpy.index_exp[:,:] + self.slice]
        else:
            return self.h5file.get_iter_group(n_iter)[self.dsname][:,:]
        
    def get_segment_data(self, n_iter, seg_id):
        if self.slice:
            return self.h5file.get_iter_group(n_iter)[numpy.index_exp[seg_id,:] + self.slice]
        else:
            return self.h5file.get_iter_group(n_iter)[seg_id,:]
        
class FnDSSpec(DSSpec):
    def __init__(self, h5file, fn):
        self.h5file = h5file
        self.fn = fn
        
    def get_iter_data(self, n_iter):
        return self.fn(n_iter, self.h5file.get_iter_group(n_iter))
        
class MultiDSSpec(DSSpec):
    def __init__(self, dsspecs):
        self.dsspecs = dsspecs
    
    def get_iter_data(self, n_iter):
        datasets = [dsspec.get_iter_data(n_iter) for dsspec in self.dsspecs]
          
        ncols = 0 
        nsegs = None
        npts = None
        for iset, dset in enumerate(datasets):
            if nsegs is None:
                nsegs = dset.shape[0]
            elif dset.shape[0] != nsegs:
                raise TypeError('dataset {} has incorrect first dimension (number of segments)'.format(self.dsspecs[iset]))
            if npts is None:
                npts = dset.shape[1]
            elif dset.shape[1] != npts:
                raise TypeError('dataset {} has incorrect second dimension (number of time points)'.format(self.dsspecs[iset]))
            
            if dset.ndim < 2:
                # scalar per segment or scalar per iteration
                raise TypeError('dataset {} has too few dimensions'.format(self.dsspecs[iset]))
            elif dset.ndim > 3:
                # array per timepoint
                raise TypeError('dataset {} has too many dimensions'.format(self.dsspecs[iset]))
            elif dset.ndim == 2:
                # scalar per timepoint
                ncols += 1
            else:
                # vector per timepoint
                ncols += dset.shape[-1]
        
        output_dtype = numpy.result_type(*[ds.dtype for ds in datasets])
        output_array = numpy.empty((nsegs, npts, ncols), dtype=output_dtype)
        
        ocol = 0
        for iset, dset in enumerate(datasets):
            if dset.ndim == 2:
                output_array[:,:,ocol] = dset[...]
                ocol += 1
            elif dset.ndim == 3:
                output_array[:,:,ocol:(ocol+dset.shape[-1])] = dset[...]
                ocol += dset.shape[-1]
        
        return output_array
    
        
    
class HistHelper:
    '''Create histograms from WEST data sets. This class coordinates reading
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
         
    '''
        
    def __init__(self, dsspec, iter_start, iter_stop, data_reader):
        self.dsspec = dsspec
        self.data_reader = data_reader
        self.iter_start = iter_start
        self.iter_stop = iter_stop
        
        self.ndim = None
        self.ntimepoints = None
        self.dset_dtype = None
        self.binbounds = None  # bin boundaries for each dimension
        self.midpoints = None  # bin midpoints for each dimension 
        self.data_range = None # data range for each dimension, as the pairs (min,max)
        self.histograms = None  # Final histogram time series
        self.avg_histogram = None # Final histogram
    
    @staticmethod    
    def parse_binspec(binspec):
        namespace = {'numpy': numpy,
                     'inf': float('inf')}
                     
        try:
            binspec_compiled = eval(binspec,namespace)
        except Exception as e:
            raise ValueError('invalid bin specification: {!r}'.format(e))
        else:
            if log.isEnabledFor(logging.DEBUG):
                log.debug('bin specs: {!r}'.format(binspec_compiled))
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
            dset = self.dsspec.get_iter_data(self.iter_start)
            self.ntimepoints = dset.shape[1]
            self.ndim = dset.shape[2]
            self.dset_dtype = dset.dtype
        
            
    def scan_data_range(self):
        '''Scan input data for range in each dimension. The number of dimensions is determined
        from the shape of the progress coordinate as of self.iter_start.'''
        
        self.scan_data_shape()
        dset_dtype = self.dset_dtype
                
        try:
            minval = numpy.finfo(dset_dtype).min
            maxval = numpy.finfo(dset_dtype).max
        except ValueError:
            minval = numpy.iinfo(dset_dtype).min
            maxval = numpy.iinfo(dset_dtype).max
        
        data_range = self.data_range = [tuple((maxval,minval)) for _i in xrange(self.ndim)]
        
        log.info('determining minimum and maximum values for {} dimensions across {} iterations'
                 .format(self.ndim,(self.iter_stop - self.iter_start)))
        for n_iter in xrange(self.iter_start, self.iter_stop):
            if log.isEnabledFor(logging.DEBUG):
                log.debug('scanning iteration {}'.format(n_iter))
            dset = self.dsspec.get_iter_data(n_iter)
            for idim in xrange(self.ndim):
                dimdata = dset[:,:,idim]
                current_min, current_max = data_range[idim]
                current_min = min(current_min, dimdata.min())
                current_max = max(current_max, dimdata.max())
                data_range[idim] = (current_min, current_max)
                del dimdata
            del dset
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
        binbounds = [numpy.require(boundset, self.dset_dtype, 'C') for boundset in self.binbounds]
        for iiter, n_iter in enumerate(xrange(self.iter_start, self.iter_stop)):
            if log.isEnabledFor(logging.DEBUG):
                log.debug('binning iteration {}'.format(n_iter))
            iter_group = self.data_reader.get_iter_group(n_iter)
            dset = self.dsspec.get_iter_data(n_iter)
            npts = dset.shape[1]
            weights = iter_group['seg_index']['weight']
            initpoint = 1 if n_iter != self.iter_start else 0
            
            dset = dset[:,initpoint:,:]
            for ipt in xrange(npts-initpoint):
                histnd(dset[:,ipt,:], binbounds, weights, out=histograms[iiter], binbound_check = False)
            
            del iter_group, weights, dset
            
            # normalize histogram
            normhistnd(histograms[iiter],binbounds)
            
        # Calculate average histogram and normalize
        # This is fast, and may be convenient for quick visualization.
        self.avg_histogram = histograms.sum(axis=0)
        normhistnd(self.avg_histogram, binbounds)
        

class WPDist(WESTTool):
    prog='w_pcpdist'
    description = '''\
Calculate probability distribution of progress coordinate values and its time evolution. 
'''
    
    def __init__(self):
        super(WPDist,self).__init__()
        
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection(self.data_reader)
        self.iter_range.include_args['iter_step'] = False
        self.binspec = None
        self.output_file = None
        self.dsspec = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
                
        #TODO: this doesn't parse right, it seems
        parser.add_argument('-b', '--bins', dest='bins', metavar='BINEXPR', default='100',
                            help='''Use BINEXPR for bins. This may be an integer, which will be used for each
                            dimension of the progress coordinate; a list of integers (formatted as [n1,n2,...])
                            which will use n1 bins for the first dimension, n2 for the second dimension, and so on;
                            or a list of lists of boundaries (formatted as [[a1, a2, ...], [b1, b2, ...], ... ]), which
                            will use [a1, a2, ...] as bin boundaries for the first dimension, [b1, b2, ...] as bin boundaries
                            for the second dimension, and so on. (Default: 100 bins in each dimension.)''')
        
        parser.add_argument('-o', '--output', dest='output', default='pdist.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        
        igroup = parser.add_argument_group('input data options').add_mutually_exclusive_group(required=True)

        igroup.add_argument('--construct-dataset',
                            help='''Use the given function (as in module.function). This function will
                            be called once per iteration as function(n_iter, iter_group) to construct
                            data for one iteration. Data returned must be indexable as
                            [seg_id][timepoint][dimension]''')
        
        igroup.add_argument('--dsspecs', nargs='+', metavar='DSSPEC',
                            help='''Construct probability distribution from DSSPEC (one per dimension).''')

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open()
        self.iter_range.process_args(args)
        
        if args.construct_dataset:
            self.dsspec = FnDSSpec(self.data_reader.we_h5file, get_object(args.construct_dataset,path=['.']))
        elif args.dsspecs:
            self.dsspec = MultiDSSpec([SingleDSSpec.from_string(dsspec, self.data_reader.we_h5file) for dsspec in args.dsspecs])
        self.binspec = args.bins        
        self.output_file = h5py.File(args.output, 'w')
        h5io.stamp_creator_data(self.output_file)
    
    def go(self):
        hh = HistHelper(self.dsspec,self.iter_range.iter_start, self.iter_range.iter_stop,self.data_reader)
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
    WPDist().main()
    
