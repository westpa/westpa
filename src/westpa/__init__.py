__all__ = [
    'Segment',
    'State',
    'Bin',
    'Source',
    'Sink',
    'Simulation',
    'TrajectoryTree',
    'TrajectoryTreeView',
    'Propagator',
    'SerialPropagator',
    'VectorizedPropagator',
    'AmberPropagator',
    'GROMACSPropagator',
    'OpenMMPropagator',
    'BinMapper',
    'RectilinearBinMapper',
    'MABBinMapper',
    'VoronoiBinMapper',
    'AdaptiveVoronoiBinMapper',
    'Resampler',
    'ResamplerBase',
    'HuberKimResampler',
    'MultinomialResampler',
    'ResidualResampler',
    'StratifiedResampler',
    'SystematicResampler',
    'SerialWorkManager',
    'ProcessWorkManager',
    'ThreadsWorkManager',
    'MPIWorkManager',
    'WESTSystem',
    'BasisState',
    'TargetState',
    '_rc',
]

import shutil

from .core.state import State
from .core.segment import Segment
from .core.protocols import Propagator, BinMapper, Resampler
from .core.propagators import SerialPropagator, VectorizedPropagator
from .core.binning import (
    Bin,
    RectilinearBinMapper,
    MABBinMapper,
    VoronoiBinMapper,
    AdaptiveVoronoiBinMapper,
)
from .core.resamplers import (
    ResamplerBase,
    HuberKimResampler,
    MultinomialResampler,
    ResidualResampler,
    StratifiedResampler,
    SystematicResampler,
)
from .core.source_sink import Source, Sink
from .core.simulation import Simulation

from .analysis import TrajectoryTree, TrajectoryTreeView
from .work_managers import SerialWorkManager, ProcessWorkManager, ThreadsWorkManager, MPIWorkManager

if shutil.which('sander'):
    from .core.propagators._amber import AmberPropagator
else:
    AmberPropagator = None

if shutil.which('gmx'):
    from .core.propagators._gromacs import GROMACSPropagator
else:
    GROMACSPropagator = None

try:
    from .core.propagators._openmm import OpenMMPropagator
except ImportError:
    OpenMMPropagator = None

from .core.states import BasisState, TargetState
from .core.systems import WESTSystem
from .core import _rc

from ._version import get_versions

rc = _rc.WESTRC()

__version__ = get_versions()['version']

del get_versions
