from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
from collections import deque
import sys
import numpy, h5py
from westtools import h5io

import westpa
from west.data_manager import weight_dtype

from westpa.kinetics import get_macrostate_rates, labeled_flux_to_rate, nested_to_flat_matrix


log = logging.getLogger('westtools.w_kinavg')


class WKinAvg(WESTTool):
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
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)

        parser.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')        
        parser.add_argument('-k', '--kinetics', default='kinetics.h5',
                            help='''Populations and transition rates are stored in KINETICS
                            (default: %(default)s).''')
        parser.add_argument('-o', '--output', dest='output', default='kinavg.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

        
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
        
    def go(self):
        nbins = self.assignments_file.attrs['nbins']
        state_labels = self.assignments_file['state_labels'][...]
        nstates = len(state_labels)
        start_iter, stop_iter = self.iter_range.iter_start, self.iter_range.iter_stop # h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter
        
        avg_flux = numpy.zeros((nstates,nstates,nbins,nbins), weight_dtype)
        avg_pops = numpy.zeros((nstates,nbins), weight_dtype)
        avg_macro_flux = numpy.zeros((nstates,nstates), weight_dtype)
        
        for iiter, n_iter in enumerate(xrange(start_iter, stop_iter)):
            if sys.stdout.isatty() and not westpa.rc.quiet_mode:
                print('\rIteration {}'.format(n_iter),end='')
                sys.stdout.flush()
            
            #assignments_iiter = h5io.get_iteration_entry(self.assignments_file, n_iter)
            kinetics_iiter = h5io.get_iteration_entry(self.kinetics_file, n_iter)
            this_flux = self.kinetics_file['labeled_bin_fluxes'][kinetics_iiter]
            this_pops = self.kinetics_file['labeled_bin_pops'][kinetics_iiter]
            this_macroflux = self.kinetics_file['trace_macro_fluxes'][kinetics_iiter]
            
            avg_flux += this_flux
            avg_pops += this_pops
            avg_macro_flux += this_macroflux            
        print()
        avg_flux /= iter_count
        avg_pops /= iter_count
        avg_macro_flux /= iter_count
        
        avg_macro_pops = avg_pops.sum(axis=1)
        
        self.output_file['avg_flux'] = avg_flux
        self.output_file['avg_flux_flat'] = nested_to_flat_matrix(avg_flux)
        
        avg_rate = labeled_flux_to_rate(avg_flux, avg_pops)
        self.output_file['avg_rate'] = avg_rate
        self.output_file['avg_rate_flat'] = nested_to_flat_matrix(avg_rate)
        ss, macro_rates = get_macrostate_rates(avg_rate, avg_pops)
        
        print(ss)
        print(state_labels)
        print(macro_rates)
        
        for istate in xrange(nstates):
            avg_macro_flux[istate,:] /= avg_macro_pops[istate]
        
        print(avg_macro_flux)
        

if __name__ == '__main__':
    WKinAvg().main()
