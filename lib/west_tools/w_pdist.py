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

from westtools import (WESTParallelTool, WESTDataReader, WESTDSSynthesizer, IterRangeSelection, 
                       ProgressIndicatorComponent)
import numpy, h5py
from fasthist import histnd, normhistnd
from westpa import h5io
from westpa.h5io import SingleIterDSSpec

log = logging.getLogger('westtools.w_pdist')

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

    data_range = [(maxval,minval) for _i in range(ndim)]

    dset = dsspec.get_iter_data(n_iter)
    for idim in range(ndim):
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
    npts = dset.shape[1]
    weights = wt_dsspec.get_iter_data(n_iter)

    dset = dset[:,initpoint:,:]
    for ipt in range(npts-initpoint):
        histnd(dset[:,ipt,:], binbounds, weights, out=iter_hist, binbound_check = False, ignore_out_of_range=ignore_out_of_range)

    del weights, dset

    # normalize histogram
    normhistnd(iter_hist,binbounds)
    return iiter, n_iter, iter_hist


class WPDist(WESTParallelTool):
    prog='w_pdist'
    description = '''\
Calculate time-resolved, multi-dimensional probability distributions of WE 
datasets.


-----------------------------------------------------------------------------
Source data
-----------------------------------------------------------------------------

Source data is provided either by a user-specified function
(--construct-dataset) or a list of "data set specifications" (--dsspecs).
If neither is provided, the progress coordinate dataset ''pcoord'' is used.

To use a custom function to extract or calculate data whose probability
distribution will be calculated, specify the function in standard Python
MODULE.FUNCTION syntax as the argument to --construct-dataset. This function
will be called as function(n_iter,iter_group), where n_iter is the iteration
whose data are being considered and iter_group is the corresponding group
in the main WEST HDF5 file (west.h5). The function must return data which can
be indexed as [segment][timepoint][dimension].

To use a list of data set specifications, specify --dsspecs and then list the
desired datasets one-by-one (space-separated in most shells). These data set
specifications are formatted as NAME[,file=FILENAME,slice=SLICE], which will
use the dataset called NAME in the HDF5 file FILENAME (defaulting to the main
WEST HDF5 file west.h5), and slice it with the Python slice expression SLICE
(as in [0:2] to select the first two elements of the first axis of the
dataset). The ``slice`` option is most useful for selecting one column (or
more) from a multi-column dataset, such as arises when using a progress
coordinate of multiple dimensions.


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
        super(WPDist,self).__init__()
        
        # Parallel processing by default (this is not actually necessary, but it is
        # informative!)
        self.wm_env.default_work_manager = self.wm_env.default_parallel_work_manager
        
        # These are used throughout
        self.progress = ProgressIndicatorComponent()
        self.data_reader = WESTDataReader()
        self.input_dssynth = WESTDSSynthesizer(default_dsname='pcoord')
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
        self.ignore_out_of_range = False
        self.compress_output = False
        
    
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
        
        parser.add_argument('-o', '--output', dest='output', default='pdist.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        parser.add_argument('-C', '--compress', action='store_true', 
                            help='''Compress histograms. May make storage of higher-dimensional histograms
                            more tractable, at the (possible extreme) expense of increased analysis time.
                            (Default: no compression.)''')
        
        parser.add_argument('--loose', dest='ignore_out_of_range', action='store_true',
                            help='''Ignore values that do not fall within bins. (Risky, as this can make buggy bin
                            boundaries appear as reasonable data. Only use if you are
                            sure of your bin boundary specification.)''')
        
        igroup = parser.add_argument_group('input dataset options').add_mutually_exclusive_group(required=False)

        igroup.add_argument('--construct-dataset',
                            help='''Use the given function (as in module.function) to extract source data.
                            This function will be called once per iteration as function(n_iter, iter_group)
                            to construct data for one iteration. Data returned must be indexable as
                            [seg_id][timepoint][dimension]''')
        
        igroup.add_argument('--dsspecs', nargs='+', metavar='DSSPEC',
                            help='''Construct probability distribution from one or more DSSPECs.''')
        
        self.progress.add_args(parser)
        
    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        self.input_dssynth.h5filename = self.data_reader.we_h5filename
        self.input_dssynth.process_args(args)
        self.dsspec = self.input_dssynth.dsspec
        
        # Carrying an open HDF5 file across a fork() seems to corrupt the entire HDF5 library
        # Open the WEST HDF5 file just long enough to process our iteration range, then close
        # and reopen in go() [which executes after the fork]
        with self.data_reader:
            self.iter_range.process_args(args)
        
        self.wt_dsspec = SingleIterDSSpec(self.data_reader.we_h5filename, 'seg_index', slice=numpy.index_exp['weight'])
        
        self.binspec = args.bins
        self.output_filename = args.output
        self.ignore_out_of_range = bool(args.ignore_out_of_range)
        self.compress_output = args.compress or False
        
    
    def go(self):
        self.data_reader.open('r')
        pi = self.progress.indicator
        pi.operation = 'Initializing'
        with pi:
            self.output_file = h5py.File(self.output_filename, 'w')
            h5io.stamp_creator_data(self.output_file)
            
            self.iter_start = self.iter_range.iter_start
            self.iter_stop = self.iter_range.iter_stop
    
            # Construct bin boundaries
            self.construct_bins(self.parse_binspec(self.binspec))
            for idim, (binbounds, midpoints) in enumerate(zip(self.binbounds, self.midpoints)):
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
        
        self.progress.indicator.new_operation('Scanning for data range', self.iter_stop-self.iter_start)
        self.scan_data_shape()
        
                
        dset_dtype = self.dset_dtype
        ndim = self.ndim
        dsspec = self.dsspec
        
        try:
            minval = numpy.finfo(dset_dtype).min
            maxval = numpy.finfo(dset_dtype).max
        except ValueError:
            minval = numpy.iinfo(dset_dtype).min
            maxval = numpy.iinfo(dset_dtype).max
        
        data_range = self.data_range = [(maxval,minval) for _i in range(self.ndim)]

        #futures = []
        #for n_iter in xrange(self.iter_start, self.iter_stop):
            #_remote_min_max(ndim, dset_dtype, n_iter, dsspec)
        #    futures.append(self.work_manager.submit(_remote_min_max, args=(ndim, dset_dtype, n_iter, dsspec)))
        
        #for future in self.work_manager.as_completed(futures):
        for future in self.work_manager.submit_as_completed(((_remote_min_max, (ndim, dset_dtype, n_iter, dsspec), {})
                                                             for n_iter in range(self.iter_start, self.iter_stop)),
                                                            self.max_queue_len):
            bounds = future.get_result(discard=True)
            for idim in range(ndim):
                current_min, current_max = data_range[idim]
                current_min = min(current_min, bounds[idim][0])
                current_max = max(current_max, bounds[idim][1])
                data_range[idim] = (current_min, current_max)
            self.progress.indicator.progress += 1

    def _construct_bins_from_scalar(self, bins):
        if self.data_range is None:
            self.scan_data_range()        

        self.binbounds = []
        self.midpoints = []        
        for idim in range(self.ndim):
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
        for idim in range(self.ndim):
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
        task_gen = ((_remote_bin_iter, (iiter, n_iter, self.dsspec, self.wt_dsspec, 1 if iiter > 0 else 0, binbounds,
                                        self.ignore_out_of_range), {}) 
                    for (iiter,n_iter) in enumerate(range(self.iter_start, self.iter_stop)))
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
    WPDist().main()
    
