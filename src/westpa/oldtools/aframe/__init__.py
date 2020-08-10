"""WEST Analyis framework -- an unholy mess of classes exploiting each other"""

from . import atool
from .atool import WESTAnalysisTool
from .base_mixin import ArgumentError, AnalysisMixin
from .binning import BinningMixin
from .data_reader import WESTDataReaderMixin, ExtDataReaderMixin, BFDataManager
from .iter_range import IterRangeMixin
from .kinetics import KineticsAnalysisMixin
from .mcbs import MCBSMixin
from .output import CommonOutputMixin
from .plotting import PlottingMixin
from .trajwalker import TrajWalker
from .transitions import TransitionAnalysisMixin, TransitionEventAccumulator, BFTransitionAnalysisMixin


__all__ = [
    'atool',
    'AnalysisMixin',
    'ArgumentError',
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
