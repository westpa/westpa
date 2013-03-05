from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTParallelTool, WESTDataReader, IterRangeSelection
import sys, math
import numpy
from westtools import h5io

import westpa
from west.data_manager import weight_dtype
import mclib
 
from westpa.kinetics import labeled_flux_to_rate, get_macrostate_rates


log = logging.getLogger('westtools.w_kinavg')

ci_dtype = numpy.dtype([('expected', numpy.float64),
                        ('ci_lbound', numpy.float64),
                        ('ci_ubound', numpy.float64)])

def _remote_get_macrostate_rates(iset, rates, pops):
    return (iset,) + get_macrostate_rates(rates,pops)

class WKinAvg(WESTParallelTool):
    prog='w_kinavg'
    description = '''\
Calculate average rates and associated errors from weighted ensemble data. Bin
assignments (usually "assignments.h5") and kinetics data (usually
"kinetics.h5") data files must have been previously generated (see 
"w_assign --help" and "w_kinetics --help" for information on generating these
files).
'''
    
    def __init__(self):
        super(WKinAvg,self).__init__()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection() 
        self.output_file = None
        self.assignments_file = None
        self.kinetics_file = None
        
        self.mcbs_alpha = None
        self.mcbs_acalpha = None
        self.mcbs_nsets = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)

        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')        
        iogroup.add_argument('-k', '--kinetics', default='kinetics.h5',
                            help='''Populations and transition rates are stored in KINETICS
                            (default: %(default)s).''')
        iogroup.add_argument('-o', '--output', dest='output', default='kinavg.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

        
        cgroup = parser.add_argument_group('confidence interval calculation options')
        cgroup.add_argument('-a', '--alpha', type=float, default=0.05, 
                             help='''Calculate a (1-ALPHA) confidence interval'
                             (default: %(default)s)''')
        cgroup.add_argument('--autocorrel-alpha', type=float, dest='acalpha', metavar='ACALPHA',
                             help='''Evaluate autocorrelation to (1-ACALPHA) significance.
                             Note that too small an ACALPHA will result in failure to detect autocorrelation
                             in a noisy flux signal. (Default: same as ALPHA.)''')
        cgroup.add_argument('-N', '--nsets', type=int,
                             help='''Use NSETS samples for bootstrapping (default: chosen based on ALPHA)''')

        
    def process_args(self, args):
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        self.kinetics_file = h5io.WESTPAH5File(args.kinetics, 'r')
        self.data_reader.process_args(args)
        self.data_reader.open('r')
        self.iter_range.process_args(args)
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        
        if not self.iter_range.check_data_iter_range_least(self.assignments_file):
            raise ValueError('assignments data do not span the requested iterations')

        if not self.iter_range.check_data_iter_range_least(self.kinetics_file):
            raise ValueError('kinetics data do not span the requested iterations')
        
        self.mcbs_alpha = args.alpha
        self.mcbs_acalpha = args.acalpha if args.acalpha else self.mcbs_alpha
        self.mcbs_nsets = args.nsets if args.nsets else mclib.get_bssize(self.mcbs_alpha)

    def kinetics_cis(self, labeled_fluxes, labeled_pops, traced_macro_fluxes):
        niters, nstates, nbins = labeled_pops.shape
        
        lbi = int(math.floor(self.mcbs_nsets*self.mcbs_alpha/2.0))
        ubi = int(math.ceil(self.mcbs_nsets*(1-self.mcbs_alpha/2.0)))

                
        ss_cidata = numpy.zeros((nstates, nbins), dtype=ci_dtype)
        mmflux_cidata = numpy.zeros((nstates,nstates), dtype=ci_dtype)
        tmflux_cidata = numpy.zeros((nstates,nstates), dtype=ci_dtype)

        # calculate means
        tmflux_cidata['expected'] = traced_macro_fluxes.mean(axis=0)        
        traj_pops = labeled_pops.mean(axis=0).sum(axis=1)
        for istate in xrange(nstates):
            if traj_pops[istate] > 0:
                tmflux_cidata['expected'][istate,:] /= traj_pops[istate]
            else:
                tmflux_cidata['expected'][istate,:] = 0

        ss_cidata['expected'], mmflux_cidata['expected'] = get_macrostate_rates(labeled_flux_to_rate(labeled_fluxes.mean(axis=0), 
                                                                                                     labeled_pops.mean(axis=0)),
                                                                                labeled_pops.mean(axis=0))
        
        # determine correlation length
        # the traced macrostate fluxes are probably a good indicator of correlation in flux data,
        # as only arrivals in each iteration are counted        
        ctimes = set()
        for istate in xrange(nstates):
            for jstate in xrange(nstates):
                if istate == jstate: continue
                ctime = mclib.mcbs_correltime(traced_macro_fluxes[:,istate,jstate], self.mcbs_acalpha, self.mcbs_nsets)
                ctimes.add(ctime)
        ctimes.discard(niters)
        if not ctimes:
            print('no correlation length determined, assuming uncorrelated fluxes')
        else:
            ctime = max(ctimes)
            print('correlation length:',ctime)
        
        stride = ctime + 1
        nslices = niters // stride
        if stride == 1:
            sliced_lfluxes = labeled_fluxes
            sliced_lpops = labeled_pops
            sliced_tmfluxes = traced_macro_fluxes
        else:
            sliced_lfluxes = numpy.empty((nslices,)+labeled_fluxes.shape[1:], labeled_fluxes.dtype)
            sliced_lpops = numpy.empty((nslices,)+labeled_pops.shape[1:], labeled_pops.dtype)
            sliced_tmfluxes = numpy.empty((nslices,)+traced_macro_fluxes.shape[1:], traced_macro_fluxes.dtype)
            for iout, istart in enumerate(xrange(0,niters-stride+1,stride)):
                sliced_lfluxes[iout] = labeled_fluxes[istart:istart+stride].mean(axis=0)
                sliced_lpops[iout] = labeled_pops[istart:istart+stride].mean(axis=0)
                sliced_tmfluxes[iout] = traced_macro_fluxes[istart:istart+stride].mean(axis=0)
                
        # We now have data blocked by correlation length; do blocked MCBS
        nsets = self.mcbs_nsets
        #synth_lfluxes = numpy.empty((nsets,)+labeled_fluxes.shape[1:], labeled_fluxes.dtype)
        #synth_lpops = numpy.empty((nsets,)+labeled_pops.shape[1:], labeled_pops.dtype)
        synth_tmfluxes = numpy.empty((nsets,)+traced_macro_fluxes.shape[1:], traced_macro_fluxes.dtype)
        #synth_lrates = numpy.empty((nsets,)+labeled_fluxes.shape[1:], labeled_fluxes.dtype)
        synth_mmfluxes = numpy.empty((nsets,)+traced_macro_fluxes.shape[1:], traced_macro_fluxes.dtype)
        synth_ss = numpy.empty((nsets,)+labeled_pops.shape[1:], labeled_pops.dtype)
        
        # Draw bootstrap data
        print('Drawing bootstrap data')
        futures = []
        for iset in xrange(self.mcbs_nsets):            
            indices = numpy.random.randint(nslices, size=(nslices,))
            synth_tmfluxes[iset] = sliced_tmfluxes[indices,...].mean(axis=0)
            synth_avgpops = sliced_lpops[indices,...].mean(axis=0)
            synth_avgflux = sliced_lfluxes[indices,...].mean(axis=0)
            synth_avgrates = labeled_flux_to_rate(synth_avgflux, synth_avgpops)
            
            # Dispatch elsewhere because this takes a while
            futures.append(self.work_manager.submit(_remote_get_macrostate_rates, args=(iset, synth_avgrates, synth_avgpops)))
            
            # Come back to local work            
            traj_pops = synth_avgpops.sum(axis=1)
            for istate in xrange(nstates):
                if traj_pops[istate] > 0:
                    synth_tmfluxes[iset,istate] /= traj_pops[istate]
                else:
                    synth_tmfluxes[iset,istate] = 0
            
            del synth_avgpops, synth_avgflux, synth_avgrates, indices, traj_pops

        nrecvd = 0        
        for future in self.work_manager.as_completed(futures):
            nrecvd += 1
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rSet {}'.format(nrecvd),end='')
                sys.stdout.flush()
            iset, ss, mmfluxes = future.get_result()
            synth_ss[iset], synth_mmfluxes[iset] = ss, mmfluxes
        print()
        
        # Calculate the (1-alpha) CI for everything
        for istate in xrange(nstates):
            for ibin in xrange(nbins):
                popelemdist = numpy.sort(synth_ss[:,istate,ibin])
                ss_cidata['ci_lbound'][istate,ibin] = popelemdist[lbi]
                ss_cidata['ci_ubound'][istate,ibin] = popelemdist[ubi]
                del popelemdist
        
        for istate in xrange(nstates):
            for jstate in xrange(nstates):
                trace_rateelemdist = numpy.sort(synth_tmfluxes[:,istate,jstate])
                matrix_rateelemdist = numpy.sort(synth_mmfluxes[:,istate,jstate])
                
                mmflux_cidata['ci_lbound'][istate,jstate] = matrix_rateelemdist[lbi]
                mmflux_cidata['ci_ubound'][istate,jstate] = matrix_rateelemdist[ubi]                
                
                tmflux_cidata['ci_lbound'][istate,jstate] = trace_rateelemdist[lbi]
                tmflux_cidata['ci_ubound'][istate,jstate] = trace_rateelemdist[ubi]
                
        
        return tmflux_cidata, mmflux_cidata, ss_cidata
        
    def go(self):
        nbins = self.assignments_file.attrs['nbins']
        state_labels = self.assignments_file['state_labels'][...]
        nstates = len(state_labels)
        start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter
                
        labeled_fluxes = numpy.empty((iter_count,nstates,nstates,nbins,nbins), weight_dtype)
        labeled_pops   = numpy.empty((iter_count,nstates,nbins), weight_dtype)
        traced_fluxes  = numpy.empty((iter_count,nstates,nstates), weight_dtype)

        print('reading kinetics data')
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rIteration {}'.format(n_iter),end='')
                sys.stdout.flush()
                
            kinetics_iiter = h5io.get_iteration_entry(self.kinetics_file, n_iter)
            labeled_fluxes[iiter] = self.kinetics_file['labeled_bin_fluxes'][kinetics_iiter]
            labeled_pops[iiter] = self.kinetics_file['labeled_bin_pops'][kinetics_iiter]
            traced_fluxes[iiter] = self.kinetics_file['trace_macro_fluxes'][kinetics_iiter]
        print()
        
        tmflux_cidata, mmflux_cidata, ss_cidata = self.kinetics_cis(labeled_fluxes, labeled_pops, traced_fluxes)
        
        self.output_file['traced_rates'] = tmflux_cidata
        self.output_file['matrix_rates'] = mmflux_cidata
        self.output_file['steady_state'] = ss_cidata
        self.output_file['state_labels'] = state_labels
        
        print('macrostate rates by tracing')
        maxlabel = max(map(len,state_labels))
        for istate in xrange(nstates):
            for jstate in xrange(nstates):
                if istate == jstate: continue
                print('{:{maxlabel}s} -> {:{maxlabel}s}: {:21.15e} ({:21.15e}, {:21.15e}) * tau^-1'
                      .format(state_labels[istate], state_labels[jstate],
                              tmflux_cidata['expected'][istate,jstate],
                              tmflux_cidata['ci_lbound'][istate,jstate],
                              tmflux_cidata['ci_ubound'][istate,jstate],
                              maxlabel=maxlabel))
        
        print('macrostate rates by steady-state extrapolation')
        for istate in xrange(nstates):
            for jstate in xrange(nstates):
                print('{:{maxlabel}s} -> {:{maxlabel}s}: {:21.15e} ({:21.15e}, {:21.15e}) * tau^-1'
                      .format(state_labels[istate], state_labels[jstate],
                              mmflux_cidata['expected'][istate,jstate],
                              mmflux_cidata['ci_lbound'][istate,jstate],
                              mmflux_cidata['ci_ubound'][istate,jstate],
                              maxlabel=maxlabel))
        

if __name__ == '__main__':
    WKinAvg().main()
