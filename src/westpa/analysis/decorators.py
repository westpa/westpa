import functools
import numpy as np

from westpa.analysis.core import Walker, Trace
from westpa.analysis.properties import WalkerProperty
from westpa.analysis.trajtools import Trajectory


def walker_property(fget=None, *, cache=True):
    """Transform a function for computing a walker property into a
    property attribute of the same name.

    Parameters
    ----------
    fget : callable
        Function for computing the property value. Must take a single
        :type:`Walker` object as input.
    cache : bool, default True
        Whether to cache the property value.

    Returns
    -------
    WalkerProperty
        The newly created property attribute. Equivalent to
        ``getattr(Walker, fget.__name__)``.

    """
    return WalkerProperty(fget, cache=cache)


def trajectory_segment(fget=None, *, cache=True):
    """Transform a function for getting a trajectory segment into a
    trajectory attribute of the same name.

    Parameters
    ----------
    fget : callable
        Function for getting a trajectory segment. Must take a single
        :type:`Walker` object as input and return a sequence.
    cached : bool, default True
        Whether to cache trajectory segments.

    Returns
    -------
    Trajectory
        The newly created trajectory attribute. Equivalent to
        ``getattr(Walker, fget.__name__)`` and
        ``getattr(Trace, fget.__name__)``.

    """
    if fget is None:
        return functools.partial(trajectory_segment, cache=cache)

    return Trajectory(fget.__name__, fget, cache_segments=cache)
