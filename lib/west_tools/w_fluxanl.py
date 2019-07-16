# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
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

import numpy, h5py
from scipy.signal import fftconvolve

from westtools import WESTTool, WESTDataReader, IterRangeSelection
import westpa
from west.data_manager import (weight_dtype, n_iter_dtype, vstr_dtype)
from west.we_driver import NewWeightEntry
import mclib
from westpa import h5io

fluxentry_dtype = numpy.dtype([('n_iter', n_iter_dtype),
                               ('flux', weight_dtype),
                               ('count', numpy.uint)])

target_index_dtype = numpy.dtype([('target_label', vstr_dtype),
                                  ('mean_flux', weight_dtype),
                                  ('mean_flux_ci_lb', weight_dtype),
                                  ('mean_flux_ci_ub', weight_dtype),
                                  ('mean_flux_correl_len', numpy.uintc)])

from westtools.dtypes import iter_block_ci_dtype as ci_dtype

def _extract_fluxes_fileversion_lt_7(iter_start, iter_stop, data_manager):
    '''Extract fluxes from old format, where groups for iterations where recyling
    occurs contain a 'recycling' table.'''
    
    assert data_manager.we_h5file_version < 7
    
    iter_count = iter_stop - iter_start
    target_count = data_manager.get_iter_group(iter_start)['recycling'].shape[0]
    fluxdata = numpy.zeros((iter_count,), dtype=fluxentry_dtype)
        
    if data_manager.we_h5file_version < 5:
        flux_field = 'weight' 
    else:
        flux_field = 'flux'
        
    fluxdata = {itarget: numpy.zeros((iter_count,), dtype=fluxentry_dtype)
                for itarget in range(target_count)}
    
    for iiter, n_iter in enumerate(range(iter_start, iter_stop)):
        rdata = data_manager.get_iter_group(n_iter)['recycling']
        for itarget in range(target_count):            
            fluxdata[itarget][iiter]['n_iter'] = n_iter
            fluxdata[itarget][iiter]['flux'] = rdata[itarget][flux_field]
            fluxdata[itarget][iiter]['count'] = rdata[itarget]['count']
        del rdata
        
    return fluxdata

def _extract_fluxes_fileversion_7(iter_start, iter_stop, data_manager):
    '''Extract fluxes from HDF5 file version 7, where recycling information is
    stored in the "new_weights" group of the iteration *following* recycling
    events.'''
    
    assert data_manager.we_h5file_version >= 7
    
    iter_count = iter_stop - iter_start
    iters = numpy.arange(iter_start, iter_stop, dtype=n_iter_dtype)
    
    # for each target by name, collect the iterations, fluxes, and counts
    # This is not the most foolproof way to do this, but it's good enough, and fast.
    # The most correct way to do this is tracing trajectories,
    # and warning if the boundary conditions change during the trace,
    # but that's for another tool.
    by_target = {}
    
    for iiter, n_iter in enumerate(range(iter_start, iter_stop)):
        target_states = data_manager.get_target_states(n_iter)
        try:
            new_weight_index = data_manager.get_iter_group(n_iter+1)['new_weights']['index']
        except KeyError:
            # no recycling data available
            continue
        
        for tstate in target_states:
            try:
                target_info = by_target[tstate.label]
            except KeyError:
                # State not seen before
                target_info = by_target[tstate.label] = numpy.zeros((iter_count,), dtype=fluxentry_dtype) 
                # If the target happens not to exist in an iteration (for whatever reason),
                # store a count of -1 as a sentinel
                target_info['count'][:] = -1
                target_info['n_iter'][:] = iters[:]
                
            recycled_from_tstate = ( (new_weight_index['source_type'] == NewWeightEntry.NW_SOURCE_RECYCLED)
                                    &(new_weight_index['target_state_id'] == tstate.state_id)
                                   )
            
            recycle_count = recycled_from_tstate.sum()
            target_info['count'][iiter] = recycle_count
            if recycle_count:
                # flux is in units of per tau
                target_info['flux'][iiter] = new_weight_index[recycled_from_tstate]['weight'].sum()
        
        del new_weight_index, target_states
    
    # Find the last contiguous run where target is available
    for target_label in by_target:
        fluxdata = by_target[target_label]
        by_target[target_label] = fluxdata[numpy.searchsorted(fluxdata['count'],[0])[0]:]
        
    return by_target

def extract_fluxes(iter_start=None, iter_stop=None, data_manager=None):
    '''Extract flux values from the WEST HDF5 file for iterations >= iter_start
    and < iter_stop, optionally using another data manager instance instead of the
    global one returned by ``westpa.rc.get_data_manager()``.
    
    Returns a dictionary mapping target names (if available, target index otherwise)
    to a 1-D array of type ``fluxentry_dtype``, which contains columns for iteration
    number, flux, and count.
    '''
    
    data_manager = data_manager or westpa.rc.get_data_manager()
    iter_start = iter_start or 1
    iter_stop = iter_stop or data_manager.current_iteration

    
    if data_manager.we_h5file_version < 7:
        return _extract_fluxes_fileversion_lt_7(iter_start, iter_stop, data_manager)
    else:
        return _extract_fluxes_fileversion_7(iter_start, iter_stop,data_manager)
    

class WFluxanlTool(WESTTool):
    prog='w_fluxanl'
    description = '''\
Extract fluxes into pre-defined target states from WEST data,
average, and construct confidence intervals. Monte Carlo bootstrapping
is used to account for the correlated and possibly non-Gaussian statistical
error in flux measurements.

All non-graphical output (including that to the terminal and HDF5) assumes that
the propagation/resampling period ``tau`` is equal to unity; to obtain results
in familiar units, divide all fluxes and multiply all correlation lengths by
the true value of ``tau``.
'''
    
    output_format_version = 2

    def __init__(self):
        super(WFluxanlTool,self).__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.output_h5file = None
        self.output_group = None
        self.target_groups = {}

        self.fluxdata = {}
        
        self.alpha = None
        self.autocorrel_alpha = None
        self.n_sets = None
        self.do_evol = False
        self.evol_step = 1
        
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)
        ogroup = parser.add_argument_group('output options')
        ogroup.add_argument('-o', '--output', default='fluxanl.h5',
                            help='Store intermediate data and analysis results to OUTPUT (default: %(default)s).')
        cgroup = parser.add_argument_group('calculation options')
        cgroup.add_argument('--disable-bootstrap', '-db', dest='bootstrap', action='store_const', const=False,
                             help='''Enable the use of Monte Carlo Block Bootstrapping.''')
        cgroup.add_argument('--disable-correl', '-dc', dest='correl', action='store_const', const=False,
                             help='''Disable the correlation analysis.''')
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
        # Disable the bootstrap or the correlation analysis.
        self.mcbs_enable = args.bootstrap if args.bootstrap is not None else True
        self.do_correl = args.correl if args.correl is not None else True
        self.autocorrel_alpha = args.acalpha or self.alpha
        self.n_sets = args.nsets or mclib.get_bssize(self.alpha)
        
        self.do_evol = args.do_evol
        self.evol_step = args.evol_step or 1
                
    def calc_store_flux_data(self):         
        westpa.rc.pstatus('Calculating mean flux and confidence intervals for iterations [{},{})'
                        .format(self.iter_range.iter_start, self.iter_range.iter_stop))
        
        fluxdata = extract_fluxes(self.iter_range.iter_start, self.iter_range.iter_stop, self.data_reader)
        
        # Create a group to store data in
        output_group = h5io.create_hdf5_group(self.output_h5file, 'target_flux', replace=False, creating_program=self.prog)        
        self.output_group = output_group
        output_group.attrs['version_code'] = self.output_format_version
        self.iter_range.record_data_iter_range(output_group)
        
        n_targets = len(fluxdata)
        index = numpy.empty((len(fluxdata),), dtype=target_index_dtype)
        avg_fluxdata = numpy.empty((n_targets,), dtype=ci_dtype)
        

        for itarget, (target_label, target_fluxdata) in enumerate(fluxdata.items()):
            # Create group and index entry
            index[itarget]['target_label'] = str(target_label)
            target_group = output_group.create_group('target_{}'.format(itarget))

            self.target_groups[target_label] = target_group
            
            # Store per-iteration values
            target_group['n_iter'] = target_fluxdata['n_iter']
            target_group['count'] = target_fluxdata['count']
            target_group['flux'] = target_fluxdata['flux']
            h5io.label_axes(target_group['flux'], ['n_iter'], units=['tau^-1'])
            
            
            # Calculate flux autocorrelation
            fluxes = target_fluxdata['flux']
            mean_flux = fluxes.mean()
            fmm = fluxes - mean_flux
            acorr = fftconvolve(fmm,fmm[::-1])
            acorr = acorr[len(acorr)//2:]
            acorr /= acorr[0]
            acorr_ds = target_group.create_dataset('flux_autocorrel', data=acorr)
            h5io.label_axes(acorr_ds, ['lag'], ['tau'])
            
            # Calculate overall averages and CIs
            #avg, lb_ci, ub_ci, correl_len = mclib.mcbs_ci_correl(fluxes, numpy.mean, self.alpha, self.n_sets,
            #                                                     autocorrel_alpha=self.autocorrel_alpha, subsample=numpy.mean)
            avg, lb_ci, ub_ci, sterr, correl_len = mclib.mcbs_ci_correl({'dataset': fluxes}, estimator=(lambda stride, dataset: numpy.mean(dataset)), alpha=self.alpha, n_sets=self.n_sets,
                                                                 autocorrel_alpha=self.autocorrel_alpha, subsample=numpy.mean, do_correl=self.do_correl, mcbs_enable=self.mcbs_enable )
            avg_fluxdata[itarget] = (self.iter_range.iter_start, self.iter_range.iter_stop, avg, lb_ci, ub_ci, sterr, correl_len)
            westpa.rc.pstatus('target {!r}:'.format(target_label))
            westpa.rc.pstatus('  correlation length = {} tau'.format(correl_len))
            westpa.rc.pstatus('  mean flux and CI   = {:e} ({:e},{:e}) tau^(-1)'.format(avg,lb_ci,ub_ci))
            index[itarget]['mean_flux'] = avg
            index[itarget]['mean_flux_ci_lb'] = lb_ci
            index[itarget]['mean_flux_ci_ub'] = ub_ci
            index[itarget]['mean_flux_correl_len'] = correl_len

        # Write index and summary        
        index_ds = output_group.create_dataset('index', data=index)
        index_ds.attrs['mcbs_alpha'] = self.alpha
        index_ds.attrs['mcbs_autocorrel_alpha'] = self.autocorrel_alpha
        index_ds.attrs['mcbs_n_sets'] = self.n_sets
        
        self.fluxdata = fluxdata
        self.output_h5file['avg_flux'] = avg_fluxdata
        
        
         
    def calc_evol_flux(self):
        westpa.rc.pstatus('Calculating cumulative evolution of flux confidence intervals every {} iteration(s)'
                        .format(self.evol_step))
        
        for itarget, (target_label, target_fluxdata) in enumerate(self.fluxdata.items()):
            fluxes = target_fluxdata['flux']
            target_group = self.target_groups[target_label]
            iter_start = target_group['n_iter'][0]
            iter_stop  = target_group['n_iter'][-1]
            iter_count = iter_stop - iter_start
            n_blocks = iter_count // self.evol_step
            if iter_count % self.evol_step > 0: n_blocks += 1
            
            cis = numpy.empty((n_blocks,), dtype=ci_dtype)
            
            for iblock in range(n_blocks):
                block_iter_stop = min(iter_start + (iblock+1)*self.evol_step, iter_stop)
                istop = min((iblock+1)*self.evol_step, len(target_fluxdata['flux']))
                fluxes = target_fluxdata['flux'][:istop]
                
                #avg, ci_lb, ci_ub, correl_len = mclib.mcbs_ci_correl(fluxes, numpy.mean, self.alpha, self.n_sets,
                #                                                     autocorrel_alpha = self.autocorrel_alpha,
                #                                                     subsample=numpy.mean)
                avg, ci_lb, ci_ub, sterr, correl_len = mclib.mcbs_ci_correl({'dataset': fluxes}, estimator=(lambda stride, dataset: numpy.mean(dataset)), alpha=self.alpha, n_sets=self.n_sets,
                                                                     autocorrel_alpha = self.autocorrel_alpha,
                                                                     subsample=numpy.mean, do_correl=self.do_correl, mcbs_enable=self.mcbs_enable )
                cis[iblock]['iter_start'] = iter_start
                cis[iblock]['iter_stop']  = block_iter_stop
                cis[iblock]['expected'], cis[iblock]['ci_lbound'], cis[iblock]['ci_ubound'] = avg, ci_lb, ci_ub
                cis[iblock]['corr_len'] = correl_len
                cis[iblock]['sterr'] = sterr
                
                del fluxes

            cis_ds = target_group.create_dataset('flux_evolution', data=cis)
            cis_ds.attrs['iter_step'] = self.evol_step
            cis_ds.attrs['mcbs_alpha'] = self.alpha
            cis_ds.attrs['mcbs_autocorrel_alpha'] = self.autocorrel_alpha
            cis_ds.attrs['mcbs_n_sets'] = self.n_sets

        
    def go(self):
        self.calc_store_flux_data()
        if self.do_evol:
            self.calc_evol_flux()
        

if __name__ == '__main__':
    WFluxanlTool().main()
