import enum
import inspect
import json
import math
from collections import UserDict

import numpy as np

from .state import State


class _AuxiliaryData(UserDict):

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise TypeError('keys must be strings, not ' + type(value).__name__)
        value = np.asarray(value)
        try:
            json.dumps(value.tolist())
        except TypeError:
            raise TypeError(f'values of dtype {str(value.dtype)!r} are not supported')
        super().__setitem__(key, value)


class Segment:
    """Data class for storing information about a trajectory segment.

    Attributes
    ----------
    n_iter : int or None
        Iteration number.
    seg_id : int or None
        Segment index.
    weight : float or None
        Statistical weight.
    parent_id : int or None
        Parent index.
    wtg_parent_ids : set of int
        Weight transfer graph parent indices.
    pcoord : 2-D numpy.ndarray or None
        Progress coordinate time series.
    status : Segment.Status or None
        Integer indicating the segment's propagation status, or None.
    initpoint_type : Segment.InitPoint
        Integer indicating the segment's origin.
    endpoint_type : Segment.EndPoint
        Integer indicating the segment's fate.
    walltime : float
        Wall-clock time taken for propagation (defaults to zero).
    cputime : float
        CPU time taken for propagation (defaults to zero).
    data : MutableMapping[str, numpy.ndarray]
        Auxiliary data.
    initial_state : State or None
        Initial state of the segment.
    final_state : State or None
        Final state of the segment.

    """

    class Status(enum.IntEnum):
        """Integer enum representing the propagation status of a segment."""

        UNSET = 0  #: Unset.
        PREPARED = 1  #: Prepared for propagation.
        COMPLETE = 2  #: Propagation complete.
        FAILED = 3  #: Propagation failed.

    class InitPoint(enum.IntEnum):
        """Integer enum representing the origin of a segment."""

        UNSET = 0  #: Unset.
        CONTINUES = 1  #: Continues a trajectory.
        NEWTRAJ = 2  #: Initiates a trajectory.

    class EndPoint(enum.IntEnum):
        """Integer enum representing the fate of a segment."""

        UNSET = 0  #: Unset.
        CONTINUES = 1  #: Trajectory continues.
        MERGED = 2  #: Trajectory pruned (merged away).
        RECYCLED = 3  #: Trajectory recycled.

    SEG_STATUS_UNSET = Status.UNSET
    SEG_STATUS_PREPARED = Status.PREPARED
    SEG_STATUS_COMPLETE = Status.COMPLETE
    SEG_STATUS_FAILED = Status.FAILED

    SEG_INITPOINT_UNSET = InitPoint.UNSET
    SEG_INITPOINT_CONTINUES = InitPoint.CONTINUES
    SEG_INITPOINT_NEWTRAJ = InitPoint.NEWTRAJ

    SEG_ENDPOINT_UNSET = EndPoint.UNSET
    SEG_ENDPOINT_CONTINUES = EndPoint.CONTINUES
    SEG_ENDPOINT_MERGED = EndPoint.MERGED
    SEG_ENDPOINT_RECYCLED = EndPoint.RECYCLED

    statuses = {f'SEG_STATUS_{member.name}': member.value for member in Status}
    initpoint_types = {f'SEG_INITPOINT_{member.name}': member.value for member in InitPoint}
    endpoint_types = {f'SEG_ENDPOINT_{member.name}': member.value for member in EndPoint}

    status_names = {member.value: f'SEG_STATUS_{member.name}' for member in Status}
    initpoint_type_names = {member.value: f'SEG_INITPOINT_{member.name}' for member in InitPoint}
    endpoint_type_names = {member.value: f'SEG_ENDPOINT_{member.name}' for member in EndPoint}

    # convenience functions for binning  # TODO: Remove.
    @staticmethod
    def initial_pcoord(segment):
        return segment.pcoord[0]

    @staticmethod
    def final_pcoord(segment):
        return segment.pcoord[-1]

    def __init__(
        self,
        n_iter=None,
        seg_id=None,
        weight=None,
        endpoint_type=None,
        parent_id=None,
        wtg_parent_ids=None,
        pcoord=None,
        status=None,
        walltime=None,
        cputime=None,
        data=None,
        initial_state=None,
        final_state=None,
    ):
        # NaNs appear sometimes if a WEST program is terminated unexpectedly; replace with zero
        walltime = 0.0 if walltime is None or math.isnan(walltime) else walltime
        cputime = 0.0 if cputime is None or math.isnan(cputime) else cputime

        # the int() and float() calls are required so that new-style string formatting doesn't barf
        # assuming that the respective fields are actually strings, probably after implicitly
        # calling __str__() on them.  Not sure if this is a numpy, h5py, or python problem
        self.n_iter = int(n_iter) if n_iter is not None else None
        self.seg_id = int(seg_id) if seg_id is not None else None
        self.status = Segment.Status(status) if status is not None else None
        self.parent_id = int(parent_id) if parent_id is not None else None
        self.endpoint_type = Segment.EndPoint(endpoint_type) if endpoint_type else Segment.EndPoint.UNSET

        self.weight = float(weight) if weight is not None else None
        self.wtg_parent_ids = set(wtg_parent_ids or ())

        self._pcoord = np.asarray(pcoord) if pcoord is not None else None
        self._walltime = walltime
        self._cputime = cputime
        self._data = _AuxiliaryData(data or {})
        self._initial_state = initial_state
        self._final_state = final_state

    @property
    def pcoord(self):
        return self._pcoord

    @pcoord.setter
    def pcoord(self, value):
        value = np.asarray(value)
        if not np.issubdtype(value.dtype, np.number):
            raise TypeError("scalar type of 'pcoord' must be numeric")
        if value.ndim != 2:
            raise ValueError("'pcoord' must be 2-D array")
        if len(value) < 2:
            raise ValueError("'pcoord' time series length must be at least 2")
        self._pcoord = value

    @property
    def walltime(self):
        return self._walltime

    @walltime.setter
    def walltime(self, value):
        value = float(value)
        if value <= 0:
            raise ValueError("'walltime' must be positive")
        self._walltime = value

    @property
    def cputime(self):
        return self._cputime

    @cputime.setter
    def cputime(self, value):
        value = float(value)
        if value <= 0:
            raise ValueError("'cputime' must be positive")
        self._cputime = value

    @property
    def data(self):
        return self._data

    @property
    def initial_state(self):
        return self._initial_state

    @initial_state.setter
    def initial_state(self, value):
        if not isinstance(value, State):
            raise TypeError("'initial_state' must be a State object")
        self._initial_state = value

    @property
    def final_state(self):
        return self._final_state

    @final_state.setter
    def final_state(self, value):
        if not isinstance(value, State):
            raise TypeError("'final_state' must be a State object")
        self._final_state = value

    def __repr__(self):
        return '<%s n_iter=%r, seg_id=%r, weight=%r, parent_id=%r at %s>' % (
            self.__class__.__name__,
            self.n_iter,
            self.seg_id,
            self.weight,
            self.parent_id,
            hex(id(self)),
        )

    @property
    def initpoint_type(self):
        if self.parent_id is None:
            return Segment.InitPoint.UNSET
        elif self.parent_id < 0:
            return Segment.InitPoint.NEWTRAJ
        else:
            return Segment.InitPoint.CONTINUES

    @property
    def initial_state_id(self):
        if self.parent_id < 0:
            return -(self.parent_id + 1)
        else:
            return None

    def __replace__(self, /, **changes):  # support copy.replace() in Python >=3.13.
        parameters = inspect.signature(self.__init__).parameters
        kwargs = {name: getattr(self, name) for name in parameters}
        kwargs.update(changes)
        return type(self)(**kwargs)

    def copy(self, **changes):
        return self.__replace__(**changes)

    # TODO: Remove. Use `segment.status.name` in new code.
    status_text = property((lambda s: s.status_names[s.status]))
    endpoint_type_text = property((lambda s: s.endpoint_type_names[s.endpoint_type]))
