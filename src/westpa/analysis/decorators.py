from westpa.analysis.core import Segment, Trace
from westpa.analysis.trajtools import Trajectory, TrajectoryConcatenator


def segment_property(fget=None, *, cache=True):
    """Transform a function for computing a segment property into a
    property attribute of the same name.

    Parameters
    ----------
    fget : callable
        Function for computing the property value. Must take a single
        :type:`Segment` object as input.
    cache : bool, default True
        Whether to cache segment property values.

    Returns
    -------
    :type:`property` or :type:`functools.cached_property`
        The newly created property attribute. Equivalent to
        ``getattr(Segment, func.__name__)``.

    """
    if fget is None:
        return functools.partial(segment_property, cache=cache)

    name = fget.__name__

    if hasattr(Segment, name):
        raise AttributeError(f'Segment.{name} already exists')

    if cache:
        prop = functools.cached_property(fget)
        prop.__set_name__(Segment, name)
    else:
        prop = property(fget)

    setattr(Segment, name, prop)

    return prop


def segment_trajectory(fget=None, *, cache=True):
    """Transform a function for getting a segment trajectory into a
    trajectory attribute of the same name.

    Parameters
    ----------
    fget : callable
        Function for getting a segment trajectory. Must take a single
        :type:`Segment` object as input and return a sequence.
    cached : bool, default True
        Whether to cache segment trajectories.

    Returns
    -------
    Trajectory
        The newly created trajectory attribute. Equivalent to both
        ``getattr(Segment, fget.__name__)`` and
        ``getattr(Trace, fget.__name__)``.

    """
    if fget is None:
        return functools.partial(segment_trajectory, cache=cache)

    name = fget.__name__

    if hasattr(Segment, name):
        raise AttributeError(f'Segment.{name} already exists')
    if hasattr(Trace, name):
        raise AttributeError(f'Trace.{name} already exists')

    concatenator = TrajectoryConcatenator(fget)
    descr = Trajectory(name, concatenator, cache_segments=cache)

    setattr(Segment, name, descr)
    setattr(Trace, name, descr)

    return descr
