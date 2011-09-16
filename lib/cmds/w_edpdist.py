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
        
    def calc_cdfs(self):
        transitions_ds = self.get_transitions_ds()
        #transitions_iter = transitions_ds['block']
        transitions_ibin = transitions_ds['initial_bin']
        transitions_fbin = transitions_ds['final_bin']
        for ibin in self.analysis_initial_bins:
            for fbin in self.analysis_final_bins:
                transitions = transitions_ds[(transitions_ibin == ibin) & (transitions_fbin == fbin)]
                wemd.rc.pstatus('  {:d}->{:d}: {:d} transitions'.format(ibin,fbin,len(transitions)),end='')
                if len(transitions):
                    edf = EDF(transitions['duration'], weights=transitions['final_weight'])
                    edf_array = edf.as_array()
                    edf_array[:,0] *= self.dt
                    wemd.rc.pstatus('; max duration = {:g}'.format(float(edf_array[-1,0])), end='')
                    numpy.savetxt('ed_{:d}_{:d}.txt'.format(ibin,fbin), edf_array)
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
wedp.calc_cdfs()
