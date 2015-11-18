from __future__ import print_function, division; __metaclass__ = type
import logging
from itertools import izip
from westtools import (WESTParallelTool, WESTDataReader, WESTDSSynthesizer, IterRangeSelection, 
                       ProgressIndicatorComponent)
import numpy, h5py
from fasthist import histnd, normhistnd
from westpa import h5io
from westpa.h5io import SingleIterDSSpec

log = logging.getLogger('westtools.w_pdist')


class DurationDataset:
    '''A facade for the 'dsspec' dataclass that incorporates the mask into get_iter_data method'''

    def __init__(self, dataset, mask, iter_start=1):
        self.dataset = dataset
        self.mask = mask
        self.dtype = dataset.dtype
        self.iter_start = iter_start

    def get_iter_data(self, n_iter):
        try:
            assert n_iter >= self.iter_start
            dset = self.dataset[n_iter-1][self.mask[n_iter-self.iter_start]]
        except(AssertionError, IndexError):
            raise ValueError, "Iteration {} is not within the iteration range".format(n_iter)
        nsegs = dset.shape[0]
        if nsegs == 0:
            return None
        else:
            return dset.reshape(nsegs, 1, 1)


def isiterable(x):
    try:
        iter(x)
    except TypeError:
        return False
    else:
        return True

def _remote_min_max(ndim, dset_dtype, n_iter, dsspec):
    try:
        minval = numpy.finfo(dset_dtype).min
        maxval = numpy.finfo(dset_dtype).max
    except ValueError:
        minval = numpy.iinfo(dset_dtype).min
        maxval = numpy.iinfo(dset_dtype).max

    data_range = [(maxval,minval) for _i in xrange(ndim)]

    dset = dsspec.get_iter_data(n_iter)

    if dset == None:
        return data_range

    for idim in xrange(ndim):
        dimdata = dset[:,:,idim]
        current_min, current_max = data_range[idim]
        current_min = min(current_min, dimdata.min())
        current_max = max(current_max, dimdata.max())
        data_range[idim] = (current_min, current_max)
        del dimdata
    del dset
    return data_range

def _remote_bin_iter(iiter, n_iter, dsspec, wt_dsspec, initpoint, binbounds, ignore_out_of_range):

    iter_hist_shape = tuple(len(bounds)-1 for bounds in binbounds)
    iter_hist = numpy.zeros(iter_hist_shape, dtype=numpy.float64)

    dset = dsspec.get_iter_data(n_iter)
    if dset is None:
        return iiter, n_iter, iter_hist
    else:
        npts = dset.shape[1]
    weights = wt_dsspec.get_iter_data(n_iter)[:,0,0]

    #dset = dset[:,initpoint:,:] 
    for ipt in xrange(npts-initpoint):
        histnd(dset[:,ipt,:], binbounds, weights, out=iter_hist, binbound_check = False, ignore_out_of_range=ignore_out_of_range)

    del weights, dset

    # normalize histogram
    normhistnd(iter_hist,binbounds)
    return iiter, n_iter, iter_hist


class WEDDist(WESTParallelTool):
    prog='w_eddist'
    description = '''\
Calculate time-resolved transition-event duration distribution from kinetics results


-----------------------------------------------------------------------------
Source data
-----------------------------------------------------------------------------

Source data is collected from the results of 'w_kinetics trace' (see w_kinetics trace --help for 
more information on generating this dataset).


-----------------------------------------------------------------------------
Histogram binning
-----------------------------------------------------------------------------

By default, histograms are constructed with 100 bins in each dimension. This
can be overridden by specifying -b/--bins, which accepts a number of different
kinds of arguments:

  a single integer N
    N uniformly spaced bins will be used in each dimension.
    
  a sequence of integers N1,N2,... (comma-separated)
    N1 uniformly spaced bins will be used for the first dimension, N2 for the
    second, and so on.
    
  a list of lists [[B11, B12, B13, ...], [B21, B22, B23, ...], ...]
    The bin boundaries B11, B12, B13, ... will be used for the first dimension,
    B21, B22, B23, ... for the second dimension, and so on. These bin
    boundaries need not be uniformly spaced. These expressions will be
    evaluated with Python's ``eval`` construct, with ``numpy`` available for
    use [e.g. to specify bins using numpy.arange()].

The first two forms (integer, list of integers) will trigger a scan of all
data in each dimension in order to determine the minimum and maximum values,
which may be very expensive for large datasets. This can be avoided by
explicitly providing bin boundaries using the list-of-lists form.

Note that these bins are *NOT* at all related to the bins used to drive WE
sampling.


-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file produced (specified by -o/--output, defaulting to "pdist.h5")
may be fed to plothist to generate plots (or appropriately processed text or
HDF5 files) from this data. In short, the following datasets are created:

  ``histograms``
    Normalized histograms. The first axis corresponds to iteration, and
    remaining axes correspond to dimensions of the input dataset.
    
  ``/binbounds_0``
    Vector of bin boundaries for the first (index 0) dimension. Additional
    datasets similarly named (/binbounds_1, /binbounds_2, ...) are created
    for additional dimensions.
    
  ``/midpoints_0``
    Vector of bin midpoints for the first (index 0) dimension. Additional
    datasets similarly named are created for additional dimensions.
    
  ``n_iter``
    Vector of iteration numbers corresponding to the stored histograms (i.e.
    the first axis of the ``histograms`` dataset).


-----------------------------------------------------------------------------
Subsequent processing
-----------------------------------------------------------------------------

The output generated by this program (-o/--output, default "pdist.h5") may be
plotted by the ``plothist`` program. See ``plothist --help`` for more
information.

    
-----------------------------------------------------------------------------
Parallelization
-----------------------------------------------------------------------------

This tool supports parallelized binning, including reading of input data.
Parallel processing is the default. For simple cases (reading pre-computed
input data, modest numbers of segments), serial processing (--serial) may be
more efficient.


-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
    
'''
    
    def __init__(self):
        super(WEDDist,self).__init__()
        
        # Parallel processing by default (this is not actually necessary, but it is
        # informative!)
        self.wm_env.default_work_manager = self.wm_env.default_parallel_work_manager
        
        # These are used throughout
        self.progress = ProgressIndicatorComponent()
        self.default_kinetics_file = 'kintrace.h5'
        self.kinetics_filename = None
        self.kinetics_file = None  #Kinavg file
        self.istate = None
        self.fstate = None
        ##Duration and weight dsspecs
        self.duration_dsspec = None
        self.wt_dsspec = None
        self.binspec = None
        self.output_filename = None
        self.output_file = None
        
        # These are used during histogram generation only
        self.iter_start = None
        self.iter_stop = None
        self.ndim = None
        #self.ntimepoints = None
        self.dset_dtype = None
        self.binbounds = None  # bin boundaries for each dimension
        self.midpoints = None  # bin midpoints for each dimension 
        self.data_range = None # data range for each dimension, as the pairs (min,max)
        self.ignore_out_of_range = False
        self.compress_output = False
        
    def add_args(self, parser):
      
        parser.add_argument('-b', '--bins', dest='bins', metavar='BINEXPR', default='100',
                            help='''Use BINEXPR for bins. This may be an integer, which will be used for each
                            dimension of the progress coordinate; a list of integers (formatted as [n1,n2,...])
                            which will use n1 bins for the first dimension, n2 for the second dimension, and so on;
                            or a list of lists of boundaries (formatted as [[a1, a2, ...], [b1, b2, ...], ... ]), which
                            will use [a1, a2, ...] as bin boundaries for the first dimension, [b1, b2, ...] as bin boundaries
                            for the second dimension, and so on. (Default: 100 bins in each dimension.)''')
        
        parser.add_argument('-C', '--compress', action='store_true', 
                            help='''Compress histograms. May make storage of higher-dimensional histograms
                            more tractable, at the (possible extreme) expense of increased analysis time.
                            (Default: no compression.)''')
        
        parser.add_argument('--loose', dest='ignore_out_of_range', action='store_true',
                            help='''Ignore values that do not fall within bins. (Risky, as this can make buggy bin
                            boundaries appear as reasonable data. Only use if you are
                            sure of your bin boundary specification.)''')

        parser.add_argument('--istate', type=int, required=True, dest='istate',
                            help='''Initial state defining transition event''')

        parser.add_argument('--fstate', type=int, required=True, dest='fstate',
                            help='''Final state defining transition event''')

        itergroup = parser.add_argument_group('iteration range options')

        itergroup.add_argument('--first-iter', default=1, dest='iter_start', type=int,
                               help='''Iteration to begin analysis (default: 1)''')

        itergroup.add_argument('--last-iter', dest='iter_stop', type=int,
                               help='''Iteration to end analysis''')
        
        iogroup = parser.add_argument_group('input/output options')
        
        # self.default_kinetics_file will be picked up as a class attribute from the appropriate subclass        
        iogroup.add_argument('-k', '--kinetics', default=self.default_kinetics_file,
                            help='''Populations and transition rates (including evolution) are stored in KINETICS
                            (default: %(default)s).''')
        iogroup.add_argument('-o', '--output', dest='output', default='eddist.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        
        self.progress.add_args(parser)
        
    def process_args(self, args):
        self.progress.process_args(args)
        self.kinetics_filename = args.kinetics
        self.istate = args.istate
        self.fstate = args.fstate
        self.kinetics_file = h5io.WESTPAH5File(self.kinetics_filename, 'r')

        self.iter_start = args.iter_start
        if args.iter_stop is None:
            self.iter_stop = self.kinetics_file.attrs['iter_stop']
        else:
            self.iter_stop = args.iter_stop + 1
        
        self.binspec = args.bins
        self.output_filename = args.output
        self.ignore_out_of_range = bool(args.ignore_out_of_range)
        self.compress_output = args.compress or False
        
    def go(self):
        
        pi = self.progress.indicator
        pi.operation = 'Initializing'
        with pi:
            self.duration = self.kinetics_file['durations'][self.iter_start-1:self.iter_stop-1]

            ##Only select transition events from specified istate to fstate
            mask = (self.duration['istate'] == self.istate) & (self.duration['fstate'] == self.fstate)

            self.duration_dsspec = DurationDataset(self.kinetics_file['durations']['duration'], mask, self.iter_start)
            self.wt_dsspec = DurationDataset(self.kinetics_file['durations']['weight'], mask, self.iter_start)

            self.output_file = h5py.File(self.output_filename, 'w')
            h5io.stamp_creator_data(self.output_file)

            # Construct bin boundaries
            self.construct_bins(self.parse_binspec(self.binspec))
            for idim, (binbounds, midpoints) in enumerate(izip(self.binbounds, self.midpoints)):
                self.output_file['binbounds_{}'.format(idim)] = binbounds
                self.output_file['midpoints_{}'.format(idim)] = midpoints

            # construct histogram
            self.construct_histogram()

            # Record iteration range        
            iter_range = numpy.arange(self.iter_start, self.iter_stop, 1, dtype=(numpy.min_scalar_type(self.iter_stop)))
            self.output_file['n_iter'] = iter_range
            self.output_file['histograms'].attrs['iter_start'] = self.iter_start
            self.output_file['histograms'].attrs['iter_stop'] = self.iter_stop
            
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
            dset = self.duration_dsspec
            #self.ntimepoints = dset.shape[1]
            #self.ndim = dset.shape[2]
            self.ndim = 1
            self.dset_dtype = dset.dtype
        
    def scan_data_range(self):
        '''Scan input data for range in each dimension. The number of dimensions is determined
        from the shape of the progress coordinate as of self.iter_start.'''
        
        self.progress.indicator.new_operation('Scanning for data range', self.iter_stop-self.iter_start)
        self.scan_data_shape()
          
        dset_dtype = self.dset_dtype
        ndim = self.ndim
        dsspec = self.duration_dsspec
        
        try:
            minval = numpy.finfo(dset_dtype).min
            maxval = numpy.finfo(dset_dtype).max
        except ValueError:
            minval = numpy.iinfo(dset_dtype).min
            maxval = numpy.iinfo(dset_dtype).max
        
        data_range = self.data_range = [(maxval,minval) for _i in xrange(self.ndim)]

        #futures = []
        #for n_iter in xrange(self.iter_start, self.iter_stop):
            #_remote_min_max(ndim, dset_dtype, n_iter, dsspec)
        #    futures.append(self.work_manager.submit(_remote_min_max, args=(ndim, dset_dtype, n_iter, dsspec)))
        
        #for future in self.work_manager.as_completed(futures):
        for future in self.work_manager.submit_as_completed(((_remote_min_max, (ndim, dset_dtype, n_iter, dsspec), {})
                                                             for n_iter in xrange(self.iter_start, self.iter_stop)),
                                                            self.max_queue_len):
            bounds = future.get_result(discard=True)
            for idim in xrange(ndim):
                current_min, current_max = data_range[idim]
                current_min = min(current_min, bounds[idim][0])
                current_max = max(current_max, bounds[idim][1])
                data_range[idim] = (current_min, current_max)
            self.progress.indicator.progress += 1

    def _construct_bins_from_scalar(self, bins):
        if self.data_range is None:
            self.scan_data_range()     

        #print(self.data_range)   

        self.binbounds = []
        self.midpoints = []        
        for idim in xrange(self.ndim):
            lb, ub = self.data_range[idim]
            # Advance just beyond the upper bound of the range, so that we catch 
            # the maximum in the histogram
            ub *= 1.01
            
            #lb -= 0.01
            
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
        The time series of histogram values is stored in ``histograms``.
        Each histogram in the time series is normalized.'''
        
        self.scan_data_shape()
        
        iter_count = self.iter_stop - self.iter_start
        histograms_ds = self.output_file.create_dataset('histograms', dtype=numpy.float64,
                                                        shape=((iter_count,) + tuple(len(bounds)-1 for bounds in self.binbounds)),
                                                        compression=9 if self.compress_output else None)
        binbounds = [numpy.require(boundset, self.dset_dtype, 'C') for boundset in self.binbounds]

        
        
        self.progress.indicator.new_operation('Constructing histograms',self.iter_stop-self.iter_start)
        task_gen = ((_remote_bin_iter, (iiter, n_iter, self.duration_dsspec, self.wt_dsspec, 0, binbounds,
                                        self.ignore_out_of_range), {}) 
                    for (iiter,n_iter) in enumerate(xrange(self.iter_start, self.iter_stop)))
        #futures = set()
        #for iiter, n_iter in enumerate(xrange(self.iter_start, self.iter_stop)):
        #    initpoint = 1 if iiter > 0 else 0
        #    futures.add(self.work_manager.submit(_remote_bin_iter,
        #                                            args=(iiter, n_iter, self.dsspec, self.wt_dsspec, initpoint, binbounds)))
        
        #for future in self.work_manager.as_completed(futures):
            #future = self.work_manager.wait_any(futures)
        #for future in self.work_manager.submit_as_completed(task_gen, self.queue_size):
        log.debug('max queue length: {!r}'.format(self.max_queue_len))
        for future in self.work_manager.submit_as_completed(task_gen, self.max_queue_len):
            iiter, n_iter, iter_hist = future.get_result(discard=True)
            self.progress.indicator.progress += 1

            # store histogram
            histograms_ds[iiter] = iter_hist
            del iter_hist, future
            

if __name__ == '__main__':
    WEDDist().main()
