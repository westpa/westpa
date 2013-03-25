from __future__ import print_function, division; __metaclass__ = type
import logging
import sys, os
from itertools import izip
from westtools.tool_classes import WESTParallelTool, WESTDataReader, IterRangeSelection
import numpy, h5py
from fasthist import histnd, normhistnd
from westpa import h5io
import westpa
from westpa.extloader import get_object
from westpa.h5io import WESTPAH5File

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
    
    def __getstate__(self):
        d = dict(self.__dict__)
        if '_h5file' in d:
            d['_h5file'] = None
        return d
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        
    
class FileLinkedDSSpec(DSSpec):
    '''Provide facilities for accessing WESTPA HDF5 files, including auto-opening and the ability
    to pickle references to such files for transmission through (e.g.) the work manager, provided
    that the HDF5 file can be accessed by the same path on both the sender and receiver.'''
    def __init__(self, h5file_or_name):
        self._h5file = None
        self._h5filename = None
        
        try:
            self._h5filename = os.path.abspath(h5file_or_name.filename)
        except AttributeError:
            self._h5filename = h5file_or_name
            self._h5file = None
        else:
            self._h5file = h5file_or_name
            
    @property
    def h5file(self):
        if self._h5file is None:
            self._h5file = WESTPAH5File(self._h5filename, 'r')
        return self._h5file
    
class SingleDSSpec(FileLinkedDSSpec):
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
        
    def __init__(self, h5file_or_name, dsname, alias=None, slice=None):
        FileLinkedDSSpec.__init__(self,h5file_or_name)
        self.dsname = dsname
        self.alias = alias or dsname
        self.slice = numpy.index_exp[slice] if slice else None

    
class SingleIterDSSpec(SingleDSSpec):
    def get_iter_data(self, n_iter):
        if self.slice:
            return self.h5file.get_iter_group(n_iter)[self.dsname][numpy.index_exp[:] + self.slice]
        else:
            return self.h5file.get_iter_group(n_iter)[self.dsname][:,:]
    
class SingleSegmentDSSpec(SingleDSSpec):
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
        
class FnDSSpec(FileLinkedDSSpec):
    def __init__(self, h5file_or_name, fn):
        FileLinkedDSSpec.__init__(self,h5file_or_name)
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
    
def _remote_bin_iter(iiter, n_iter, dsspec, wt_dsspec, initpoint, binbounds):
    
    iter_hist_shape = tuple(len(bounds)-1 for bounds in binbounds)
    iter_hist = numpy.zeros(iter_hist_shape, dtype=numpy.float64)

    dset = dsspec.get_iter_data(n_iter)
    npts = dset.shape[1]
    weights = wt_dsspec.get_iter_data(n_iter)
    
    dset = dset[:,initpoint:,:]
    for ipt in xrange(npts-initpoint):
        histnd(dset[:,ipt,:], binbounds, weights, out=iter_hist, binbound_check = False)
    
    del weights, dset
    
    # normalize histogram
    normhistnd(iter_hist,binbounds)
    return iiter, n_iter, iter_hist

    
class WPDist(WESTParallelTool):
    prog='w_pdist'
    description = '''\
Calculate probability distribution of progress coordinate values and its time evolution. 
'''
    
    def __init__(self):
        super(WPDist,self).__init__()
        
        # These are used throughout
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection(self.data_reader)
        self.iter_range.include_args['iter_step'] = False
        self.binspec = None
        self.output_filename = None
        self.output_file = None
        self.dsspec = None
        self.wt_dsspec = None # dsspec for weights
        
        # These are used during histogram generation only
        self.iter_start = None
        self.iter_stop = None
        self.ndim = None
        self.ntimepoints = None
        self.dset_dtype = None
        self.binbounds = None  # bin boundaries for each dimension
        self.midpoints = None  # bin midpoints for each dimension 
        self.data_range = None # data range for each dimension, as the pairs (min,max)
        
    
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
        
        igroup = parser.add_argument_group('input data options').add_mutually_exclusive_group(required=False)

        igroup.add_argument('--construct-dataset',
                            help='''Use the given function (as in module.function). This function will
                            be called once per iteration as function(n_iter, iter_group) to construct
                            data for one iteration. Data returned must be indexable as
                            [seg_id][timepoint][dimension]''')
        
        igroup.add_argument('--dsspecs', nargs='+', metavar='DSSPEC',
                            help='''Construct probability distribution from DSSPEC (one dimension per column).''')

    def process_args(self, args):
        self.data_reader.process_args(args)
        
        # Carrying an open HDF5 file across a fork() seems to corrupt the entire HDF5 library
        # Open the WEST HDF5 file just long enough to process our iteration range, then close
        # and reopen in go() [which executes after the fork]
        self.data_reader.open('r')
        self.iter_range.process_args(args)
        self.data_reader.close()
        
        if args.construct_dataset:
            self.dsspec = FnDSSpec(self.data_reader.we_h5filename, get_object(args.construct_dataset,path=['.']))
        elif args.dsspecs:
            self.dsspec = MultiDSSpec([SingleSegmentDSSpec.from_string(dsspec, self.data_reader.we_h5filename)
                                       for dsspec in args.dsspecs])
        else:
            self.dsspec = SingleSegmentDSSpec(self.data_reader.we_h5filename, 'pcoord')
            
        self.wt_dsspec = SingleIterDSSpec(self.data_reader.we_h5filename, 'seg_index', slice=numpy.index_exp['weight'])
            
        self.binspec = args.bins
        self.output_filename = args.output
        
    
    def go(self):
        self.data_reader.open('r')
        self.output_file = h5py.File(self.output_filename, 'w')
        h5io.stamp_creator_data(self.output_file)
        
        self.iter_start = self.iter_range.iter_start
        self.iter_stop = self.iter_range.iter_stop

        # Construct bin boundaries
        self.construct_bins(self.parse_binspec(self.binspec))
        for idim, (binbounds, midpoints) in enumerate(izip(self.binbounds, self.midpoints)):
            self.output_file['binbounds_{}'.format(idim)] = binbounds
            self.output_file['midpoints_{}'.format(idim)] = midpoints

        # construct histogram
        self.construct_histogram()

        # Record iteration range        
        iter_range = self.iter_range.iter_range()
        self.output_file['n_iter'] = iter_range
        self.iter_range.record_data_iter_range(self.output_file['histograms'])
        
        self.output_file.close()

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
        
        if sys.stdout.isatty() and not westpa.rc.quiet_mode:
            print('Scanning for minimum/maximum values')
            
        for n_iter in xrange(self.iter_start, self.iter_stop):
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rIteration {}'.format(n_iter), end='')
            
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
        if sys.stdout.isatty() and not westpa.rc.quiet_mode:
            print('')
            
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
        histograms_ds = self.output_file.create_dataset('histograms', dtype=numpy.float64,
                                                        shape=((iter_count,) + tuple(len(bounds)-1 for bounds in self.binbounds)))
        binbounds = [numpy.require(boundset, self.dset_dtype, 'C') for boundset in self.binbounds]
        
        if sys.stdout.isatty() and not westpa.rc.quiet_mode:
            print('Creating histograms')
        
        futures = []
        for iiter, n_iter in enumerate(xrange(self.iter_start, self.iter_stop)):                
            initpoint = 1 if iiter > 0 else 0
            futures.append(self.work_manager.submit(_remote_bin_iter,
                                                    args=(iiter, n_iter, self.dsspec, self.wt_dsspec, initpoint, binbounds)))
        
        n_received = 0
        for future in self.work_manager.as_completed(futures):
            iiter, n_iter, iter_hist = future.get_result(discard=True)
            n_received += 1

            # store histogram
            histograms_ds[iiter] = iter_hist
            del iter_hist, future

            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rFinished {} of {} iterations'.format(n_received,iter_count), end='')            
        if sys.stdout.isatty() and not westpa.rc.quiet_mode:
            print('')


if __name__ == '__main__':
    WPDist().main()
    
