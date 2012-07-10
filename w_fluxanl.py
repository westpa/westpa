from __future__ import print_function, division; __metaclass__ = type
import sys
from wt2.tool_classes import WEMDTool, HDF5Storage, WEMDDataReader, IterRangeSelection
from itertools import imap
import numpy, h5py, operator, functools
import scipy.signal
from scipy.signal import fftconvolve

import wemd
from wemd.data_manager import (weight_dtype, n_iter_dtype)

import mclib

ci_dtype = numpy.dtype([('mean', numpy.float64),
                        ('ci_lb', numpy.float64),
                        ('ci_ub', numpy.float64),
                        ('correl_len', numpy.uintc)])

iter_range_dtype = numpy.dtype([('iter_start', n_iter_dtype),
                                ('iter_stop', n_iter_dtype)])


def extract_fluxes(iter_start=None, iter_stop=None, data_manager=None):
    '''Extract flux values from the WEMD HDF5 file for iterations >=iter_start
    and <iter_stop, optionally using another data manager instance instead of the
    global one returned by ``wemd.rc.get_data_manager()``.  Returns a triplet
    ``(iters, fluxes, counts)`` where ``iters`` is an array of the iterations
    considered, ``fluxes`` is an array of flux values, indexed as
    ``fluxes[n_iter][itarget]``, and ``counts`` is an array of recycling 
    counts, indexed as ``counts[n_iter][itarget]``.'''
    
    data_manager = data_manager or wemd.rc.get_data_manager()
    iter_start = iter_start or 1
    iter_stop = iter_stop or data_manager.current_iteration
    iter_count = iter_stop - iter_start
    target_count = data_manager.get_iter_group(iter_start)['recycling'].shape[0]
    
    iters = numpy.arange(iter_start, iter_stop, dtype=n_iter_dtype)
    fluxes = numpy.zeros((iter_count, target_count), weight_dtype)
    counts = numpy.zeros((iter_count, target_count), numpy.uint)
    
    if data_manager.we_h5file_version < 5:
        flux_field = 'weight' 
    else:
        flux_field = 'flux'
    
    for iiter, n_iter in enumerate(xrange(iter_start, iter_stop)):
        rdata = data_manager.get_iter_group(n_iter)['recycling']
        
        fluxes[iiter,:] = rdata[flux_field]
        counts[iiter,:] = rdata['count']
        
    return (iters,fluxes,counts)

class WFluxanlTool(WEMDTool):
    prog='w_fluxanl'
    description = '''\
Extract fluxes into pre-defined target states from WEMD data,
average, and construct confidence intervals. Monte Carlo bootstrapping
is used to account for the correlated and possibly non-Gaussian statistical
error in flux measurements.

All non-graphical output (including that to the terminal and HDF5) assumes that
the propagation/resampling period ``tau`` is equal to unity; to obtain results
in familiar units, divide all fluxes and multiply all correlation lengths by
the true value of ``tau``.
'''

    def __init__(self):
        super(WFluxanlTool,self).__init__()
        self.data_reader = WEMDDataReader()
        self.iter_range = IterRangeSelection()
        self.output_h5file = None
        
        self.alpha = None
        self.autocorrel_alpha = None
        self.n_sets = None
        
        self.iters = None
        self.fluxes = None
        self.counts = None
        self.n_targets = None
        
        self.do_evol = False
        self.evol_step = 1
        
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
        ogroup = parser.add_argument_group('output options')
        ogroup.add_argument('-o', '--output', default='fluxanl.h5',
                            help='Store intermediate data and analysis results to OUTPUT (default: %(default)s).')
        cgroup = parser.add_argument_group('calculation options')
        cgroup.add_argument('-a', '--alpha', type=float, default=0.05, 
                             help='''Calculate a (1-ALPHA) confidence interval on the average flux'
                             (default: %(default)s)''')
        cgroup.add_argument('--autocorrel-alpha', type=float, dest='acalpha', metavar='ACALPHA',
                             help='''Evaluate autocorrelation of flux to (1-ACALPHA) significance.
                             Note that too small an ACALPHA will result in failure to detect autocorrelation
                             in a noisy flux signal. (Default: same as ALPHA.)''')
        cgroup.add_argument('-N', '--nsets', type=int,
                             help='''Use NSETS samples for bootstrapping (default: chosen based on ALPHA)''')
        cgroup.add_argument('--evol', action='store_true', dest='do_evol',
                            help='''Calculate time evolution of flux confidence intervals (expensive).''')
        cgroup.add_argument('--evol-step', type=int, default=1, metavar='ESTEP',
                            help='''Calculate time evolution of flux confidence intervals every ESTEP
                            iterations (default: %(default)s)''')
        
        
    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open()
        self.iter_range.data_manager = self.data_reader
        self.iter_range.process_args(args)
        
        self.output_h5file = h5py.File(args.output, 'w')
        
        self.alpha = args.alpha
        self.autocorrel_alpha = args.acalpha or self.alpha
        self.n_sets = args.nsets or mclib.get_bssize(self.alpha)
        
        self.do_evol = args.do_evol
        self.evol_step = args.evol_step or 1
                
    def calc_store_flux_data(self):         
        wemd.rc.pstatus('Extracting fluxes and transition counts for iterations [{},{})'
                        .format(self.iter_range.iter_start, self.iter_range.iter_stop))
        
        iters, fluxes, counts = extract_fluxes(self.iter_range.iter_start, self.iter_range.iter_stop, self.data_reader)
        self.output_h5file['iterations'] = iters
        self.output_h5file['fluxes'] = fluxes
        self.output_h5file['fluxes'].attrs['description'] = 'instantaneous flux'
        self.output_h5file['counts'] = counts
        self.output_h5file['counts'].attrs['description'] = 'instantaneous transition counts'
        
        # Stamp data sets with axis labels
        for h5object in (self.output_h5file['fluxes'], self.output_h5file['counts']):
            h5object.attrs['axis_labels'] = numpy.array(['n_iter','target_index'])

        # Stamp data sets with range over which this analysis was conducted
        for h5object in (self.output_h5file, self.output_h5file['iterations'],
                         self.output_h5file['fluxes'], self.output_h5file['counts']):
            self.iter_range.record_data_iter_range(h5object)
            
        self.iters = iters
        self.fluxes = fluxes
        self.counts = counts
        self.n_targets = self.fluxes.shape[1]
        
    def calc_store_flux_autocorrel(self):
        fmm = self.fluxes - self.fluxes.mean(axis=0)
        acorr = numpy.empty((self.fluxes.shape[0], self.fluxes.shape[1]), numpy.float64)
        
        for target in xrange(self.n_targets):
            target_acorr = fftconvolve(fmm[:,target], fmm[::-1,target])
            target_acorr = target_acorr[len(target_acorr)//2:]
            target_acorr /= target_acorr[0]
            acorr[:,target] = target_acorr
            del target_acorr
            
        self.output_h5file['autocorrel'] = acorr
        h5ds = self.output_h5file['autocorrel']
        h5ds.attrs['description'] = 'flux autocorrelation'
        h5ds.attrs['axis_labels'] = numpy.array(['lag','target_index'])        
        self.iter_range.record_data_iter_range(h5ds)
        
    def calc_overall_avg_flux(self):
        wemd.rc.pstatus('Calculating alpha={} confidence interval on mean flux for {} target states'
                        .format(self.alpha, self.n_targets))
        
        cis = numpy.empty((self.n_targets,), ci_dtype)
        
        for target in xrange(self.n_targets):
            cis[target]   = avg, lb_ci, ub_ci, correl_len \
                          = mclib.mcbs_ci_correl(self.fluxes[:,target], numpy.mean, self.alpha, self.n_sets,
                                                 autocorrel_alpha=self.autocorrel_alpha, subsample=numpy.mean)
            wemd.rc.pstatus('  target {}:'.format(target))
            wemd.rc.pstatus('    correlation length = {} tau'.format(correl_len))
            wemd.rc.pstatus('    mean flux and CI   = {} ({},{}) tau^(-1)'.format(avg, lb_ci, ub_ci))
            
        h5ds = self.output_h5file.create_dataset('overall_average', data=cis)
        self.iter_range.record_data_iter_range(h5ds)
        h5ds.attrs['description'] = 'flux mean and confidence interval over entire data range'
        h5ds.attrs['axis_labels'] = numpy.array(['target_index'])
        h5ds.attrs['mcbs_alpha'] = self.alpha
        h5ds.attrs['mcbs_autocorrel_alpha'] = self.autocorrel_alpha
        h5ds.attrs['mcbs_n_sets'] = self.n_sets
         
    def calc_evol_flux(self):
        wemd.rc.pstatus('Calculating cumulative evolution of flux confidence intervals every {} iteration(s)'
                        .format(self.evol_step))
        
        
        iter_start = self.iter_range.iter_start
        n_blocks = self.iter_range.iter_count // self.evol_step
        if self.iter_range.iter_count % self.evol_step > 0:
            n_blocks += 1
        
        iters = numpy.empty((n_blocks,), dtype=iter_range_dtype)
        cis = numpy.empty((n_blocks,self.n_targets), dtype=ci_dtype)
        
        for iblock in xrange(0,n_blocks):
            iter_stop = min(iter_start + (iblock+1)*self.evol_step, self.iter_range.iter_stop)
            istop = min((iblock+1)*self.evol_step, len(self.fluxes))
            iters[iblock]['iter_start'] = iter_start
            iters[iblock]['iter_stop']  = iter_stop
            
            for target in xrange(self.n_targets):
                fluxes = self.fluxes[:istop,target]
                cis[iblock,target] = mclib.mcbs_ci_correl(fluxes, numpy.mean, self.alpha, self.n_sets,
                                                          autocorrel_alpha=self.autocorrel_alpha, subsample=numpy.mean)
                del fluxes
        
        self.output_h5file.create_dataset('flux_evol_iterations', data=iters)
        evol_cis_ds = self.output_h5file.create_dataset('flux_evol', data=cis)
        self.iter_range.record_data_iter_range(evol_cis_ds)
        evol_cis_ds.attrs['iter_step'] = self.evol_step
        evol_cis_ds.attrs['description'] = 'time-resolved cumulative flux mean and confidence intervals'
        evol_cis_ds.attrs['axis_labels'] = numpy.array(['iteration_block', 'target_index'])
        evol_cis_ds.attrs['mcbs_alpha'] = self.alpha
        evol_cis_ds.attrs['mcbs_autocorrel_alpha'] = self.autocorrel_alpha
        evol_cis_ds.attrs['mcbs_n_sets'] = self.n_sets
        
    def go(self):
        self.calc_store_flux_data()
        self.calc_store_flux_autocorrel()
        self.calc_overall_avg_flux()
        if self.do_evol:
            self.calc_evol_flux()
        

if __name__ == '__main__':
    WFluxanlTool().main()