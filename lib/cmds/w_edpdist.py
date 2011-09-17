from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse
import numpy
import wemd
from wemdtools.stats.edfs import EDF

import logging
log = logging.getLogger('w_edpdist')

from wemdtools.aframe import (WEMDAnalysisTool,WEMDDataReaderMixin,IterRangeMixin,TransitionAnalysisMixin,BinningMixin,
                              KineticsAnalysisMixin)
                              
class WEDPDist(KineticsAnalysisMixin,TransitionAnalysisMixin,BinningMixin,IterRangeMixin,WEMDDataReaderMixin,WEMDAnalysisTool):
    def __init__(self):
        super(WEDPDist,self).__init__()
        self.edp_group = None
        self.edfs_group = None
        
    def require_edf_group(self):
        self.edfs_group = self.edp_group.require_group('edfs')
        for (k,v) in self.get_transitions_ds().attrs.iteritems():
            self.edfs_group.attrs[k] = v
        
    def store_edf(self, ibin, fbin, first_iter, last_iter, edf):
        '''Store the given EDF, which is for ibin->fbin transitions, in the HDF5 file.'''
        ibin_group = self.edfs_group.require_group(str(ibin))
        ibin_group.attrs['initial_bin'] = ibin
        fbin_group = ibin_group.require_group(str(fbin))
        fbin_group.attrs['final_bin'] = fbin
        fbin_group['edf_iters_{:d}_{:d}'.format(first_iter, last_iter)] = edf.as_array()
        
    def calc_cdfs(self):
        wemd.rc.pstatus('Calculating event duration CDFs...')
        transitions_ds = self.get_transitions_ds()
        transitions_niter = transitions_ds['n_iter']
        transitions_ibin = transitions_ds['initial_bin']
        transitions_fbin = transitions_ds['final_bin']
        for ibin in self.analysis_initial_bins:
            for fbin in self.analysis_final_bins:
                tdisc = (transitions_ibin == ibin) & (transitions_fbin == fbin)
                tdisc &= (transitions_niter >= self.first_iter) & (transitions_niter <= self.last_iter)
                transitions = transitions_ds[tdisc]
                wemd.rc.pstatus('  {:d}->{:d}: {:d} transitions'.format(ibin,fbin,len(transitions)),end='')
                if len(transitions):
                    durations = transitions['duration'] * self.dt
                    weights = transitions['final_weight']
                    edf = EDF(durations, weights)
                    edf_array = edf.as_array()
                    wemd.rc.pstatus('; N_e = {:d}, mean = {:g}, stdev = {:g}, median = {:g}, 95th percentile = {:g}, max = {:g}'
                                    .format(len(edf),
                                            float(edf.mean()),
                                            float(edf.std()),
                                            float(edf.quantile(0.5)),
                                            float(edf.quantile(0.95)),
                                            float(edf_array[-1,0])), 
                                    end='')
                    self.store_edf(ibin, fbin, self.first_iter, self.last_iter, edf)
                    #numpy.savetxt('ed_{:d}_{:d}.txt'.format(ibin,fbin), edf_array)
                    del edf
                del transitions
                wemd.rc.pstatus()


        
wedp = WEDPDist()

parser = argparse.ArgumentParser('w_edpdist', description='''\
Calculate the transition event duration probability distribution as a function of number of iterations,
and evaluate WEMD simulation convergence. 
''')
wemd.rc.add_args(parser)
wedp.add_args(parser)

args = parser.parse_args()

wemd.rc.process_args(args, config_required=False)
wedp.process_args(args)
wedp.dt = args.dt

# Preliminary checks and required calculation steps
wedp.check_iter_range()
wedp.check_bin_selection()
wedp.open_analysis_backing()
wedp.edp_group = wedp.require_analysis_group('w_edpdist', replace=True)
wedp.require_bin_assignments()
wedp.require_transitions()
wedp.require_edf_group()
wedp.calc_cdfs()
