# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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

"""WEST Analyis framework -- an unholy mess of classes exploiting each other"""


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
