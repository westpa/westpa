"""WEST Analyis framework -- an unholy mess of classes exploiting each other"""

from . import atool
from .atool import WESTAnalysisTool
from .iter_range import IterRangeMixin
from .data_reader import WESTDataReaderMixin, ExtDataReaderMixin, BFDataManager
from .binning import BinningMixin
from .mcbs import MCBSMixin
from .trajwalker import TrajWalker
from .transitions import TransitionAnalysisMixin, TransitionEventAccumulator, BFTransitionAnalysisMixin
from .kinetics import KineticsAnalysisMixin
from .output import CommonOutputMixin
from .plotting import PlottingMixin


class ArgumentError(RuntimeError):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class AnalysisMixin:
    def __init__(self):
        super().__init__()

    def add_args(self, parser, upcall=True):
        if upcall:
            try:
                upfunc = super().add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)

    def process_args(self, args, upcall=True):
        if upcall:
            try:
                upfunc = super().process_args
            except AttributeError:
                pass
            else:
                upfunc(args)


__all__ = [
    'atool',
    'WESTAnalysisTool',
    'IterRangeMixin',
    'WESTDataReaderMixin',
    'ExtDataReaderMixin',
    'BFDataManager',
    'BinningMixin',
    'MCBSMixin',
    'TrajWalker',
    'TransitionAnalysisMixin',
    'TransitionEventAccumulator',
    'BFTransitionAnalysisMixin',
    'KineticsAnalysisMixin',
    'CommonOutputMixin',
    'PlottingMixin',
]
