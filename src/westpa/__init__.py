from ._version import get_versions

from .core.segment import Segment
from .core.systems import WESTSystem
from .core.states import BasisState, TargetState
from .core import _rc

__all__ = ['Segment', 'WESTSystem', 'BasisState', 'TargetState', '_rc']

rc = _rc.WESTRC()

__version__ = get_versions()["version"]

del get_versions
