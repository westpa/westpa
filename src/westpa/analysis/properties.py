import functools
import operator

from westpa.analysis.core import Segment, Trace


def segment_property(fget=None, *, cached=True):
    def register_segment_property(fget):
        name = fget.__name__

        if hasattr(Segment, name):
            raise AttributeError(f'Segment.{name} already exists')

        if cached:
            prop = functools.cached_property(fget)
            prop.__set_name__(Segment, name)
        else:
            prop = property(fget)

        setattr(Segment, name, prop)

        return fget

    if fget is None:
        return register_segment_property
    return register_segment_property(fget)


class SegmentTrajectory:

    def __init__(self, fget, name, cached=True):
        self.fget = fget
        self.name = name
        self.cached = cached

    @property
    def private_name(self):
        return '_' + self.name

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self

        if hasattr(obj, self.private_name):
            value = getattr(obj, self.private_name)
        else:
            value = self.fget(obj)

        if not hasattr(value, '__getitem__'):
            msg = f"'{type(value).__name__}' object can't be concatenated"
            raise TypeError(msg)

        if self.cached:
            setattr(obj, self.private_name, value)

        return value

    def __set__(self, obj, value):
        raise AttributeError("can't set attribute")


class TraceTrajectory:

    def __init__(self, name):
        self.name = name

    def get_segment_traj(self, segment):
        return getattr(segment, self.name)

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        segment_trajs = map(self.get_segment_traj, obj.segments)
        return functools.reduce(operator.concat, segment_trajs)

    def __set__(self, obj, value):
        raise AttributeError("can't set attribute")


def segment_trajectory(fget=None, *, cached=True, concat_operator=None):
    def register_segment_trajectory(fget):
        name = fget.__name__

        if hasattr(Segment, name):
            raise AttributeError(f'Segment.{name} already exists')
        if hasattr(Trace, name):
            raise AttributeError(f'Trace.{name} already exists')

        setattr(Segment, name, SegmentTrajectory(fget, name, cached=cached))
        setattr(Trace, name, TraceTrajectory(name))

        return fget

    if fget is None:
        return register_segment_trajectory
    return register_segment_trajectory(fget)
