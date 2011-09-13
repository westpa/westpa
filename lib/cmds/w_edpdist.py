from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse
import numpy
import wemd

import logging
log = logging.getLogger('w_edpdist')

from wemdtools.aframe import (WEMDAnalysisTool,DataReaderMixin,IterRangeMixin,TransitionAnalysisMixin,BinningMixin,
                              KineticsAnalysisMixin)
                              
class WEDPDist(KineticsAnalysisMixin,TransitionAnalysisMixin,BinningMixin,IterRangeMixin,DataReaderMixin,WEMDAnalysisTool):
    def __init__(self):
        super(WEDPDist,self).__init__()
        
    def calc_cdfs(self):
        transitions_ds = self.get_transitions_ds()
        transition_iters = transitions_ds['block']
        for (block_first_iter,block_past_last) in self.iter_block_iter():
            wemd.rc.pstatus('\r  Iterations [{:d},{:d}): '.format(block_first_iter,block_past_last), end='')
            wemd.rc.pflush()
            transitions = transitions_ds[(transition_iters >= block_first_iter) & (transition_iters < block_past_last)]
            transitions = transitions[numpy.in1d(transitions['initial_bin'], self.analysis_initial_bins)]
            transitions = transitions[numpy.in1d(transitions['final_bin'], self.analysis_final_bins)]
            
            wemd.rc.pflush()
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
