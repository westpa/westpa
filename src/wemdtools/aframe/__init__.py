"""WEMD Analyis framework"""

from __future__ import division, print_function; __metaclass__ = type

class ArgumentError(RuntimeError):
    def __init__(self, *args, **kwargs):
        super(ArgumentError,self).__init__(*args,**kwargs)

class AnalysisMixin:
    def __init__(self):
        super(AnalysisMixin,self).__init__()
        
    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upfunc = super(AnalysisMixin,self).add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)
    
    def process_args(self, args, upcall = True):
        if upcall:
            try:
                upfunc = super(AnalysisMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)
    
import atool
from atool import WEMDAnalysisTool
from iter_range import IterRangeMixin
from data_reader import WEMDDataReaderMixin
from binning import BinningMixin
from mcbs import MCBSMixin
from trajwalker import TrajWalker
from transitions import TransitionAnalysisMixin, TransitionEventAccumulator
from kinetics import KineticsAnalysisMixin
from output import CommonOutputMixin
from plotting import PlottingMixin
