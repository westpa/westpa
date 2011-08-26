from __future__ import division, print_function
import os, sys, argparse, numpy
from math import ceil, log10

import logging
log = logging.getLogger('w_binprobs')

import wemd, wemdtools

from wemdtools.aframe import WEMDAnalysisTool, BinningMixin, DataReaderMixin, IterRangeMixin

# Upcalls move left to right
# If some complex dependencies exist, one can always override the process_common_args function and
# call parent classes' process_common_args() manually in an order that makes sense.
class WBinprobs(BinningMixin, DataReaderMixin, IterRangeMixin, WEMDAnalysisTool):
    pass

wbp = WBinprobs()

parser = argparse.ArgumentParser('w_binprobs')
wemd.rc.add_common_args(parser)
wbp.add_common_args(parser)

args = parser.parse_args()

wemd.rc.process_common_args(args, config_required=False)
wbp.process_common_args(args)
wbp.check_iter_range()
wbp.open_analysis_backing()
wbp.require_bin_assignments()
