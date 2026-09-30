import enum
import inspect
from collections import UserDict
from collections.abc import MutableSet

import numpy as np

from .state import State


# Used to validate the 'data' attribute.
# Ensures that keys are strings and values are non-object arrays.
class _AuxiliaryData(UserDict):

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise TypeError('keys must be strings, not ' + type(value).__name__)
        value = np.asarray(value)
        if value.dtype == object:
            raise TypeError('object arrays are not supported')
        super().__setitem__(key, value)


# Used to validate the 'wtg_parent_ids' attribute.
class _IntegerSet(MutableSet):

    def __init__(self, iterable=None):
        self._data = set()
        for elem in iterable or ():
            self.add(elem)

    def __repr__(self):
        return repr(self._data)

    def __contains__(self, item):
        return item in self._data

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def add(self, elem):
        if not isinstance(elem, int):
            raise TypeError(f"elements must be integers, not {type(elem).__name__}")
        self._data.add(elem)

    def discard(self, elem):
        self._data.discard(elem)


class Segment:
    """Represents a trajectory segment.

    Parameters
    ----------
    n_iter : int, optional
        Iteration number.
    seg_id : int, optional
        Segment index.
    weight : float, optional
        Statistical weight.
    endpoint_type : Segment.EndPoint, optional
        Integer indicating whether the trajectory continues to the next
        iteration, gets recycled, or gets merged away. Defaults to
        ``<EndPoint.UNSET: 0>``.
    parent_id : int, optional
        Index of the segment's parent. A negative value indicates a new trajectory.
    wtg_parent_ids : iterable of int, optional
        Indices of the walkers from which the segment inherited weight.
    pcoord : 2-D array_like, optional
        Progress coordinate time series. The time series must include at least
        the initial and final points (the array length must be at least two).
    status : Segment.Status, optional
        Integer indicating the propagation status of the segment. Defaults to
        ``<Status.UNSET: 0>``.
    walltime : float, optional
        Wall-clock time (in seconds) taken for propagation. Defaults to zero.
    cputime : float, optional
        CPU time (in seconds) taken for propagation. Defaults to zero.
    data : Mapping[str, array_like], optional
        Auxiliary data. Object arrays are not supported.
    initial_state : State, optional
        Initial state of the segment.
    final_state : State, optional
        Final state of the segment.

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
    wtg_parent_ids : MutableSet[int]
        Weight transfer graph parent indices.
    pcoord : numpy.ndarray or None
        Progress coordinate time series.
    status : Segment.Status
        Propagation status.
    initpoint_type : Segment.InitPoint
        Whether the segment continues a trajectory or starts a new trajectory.
    endpoint_type : Segment.EndPoint
        Whether the trajectory continues, gets recycled, or gets merged away.
    walltime : float
        Wall-clock time taken for propagation.
    cputime : float
        CPU time taken for propagation.
    data : MutableMapping[str, numpy.ndarray]
        Auxiliary data.
    initial_state : State or None
        Initial state.
    final_state : State or None
        Final state.

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
        NEWTRAJ = 2  #: Starts a new trajectory.

    class EndPoint(enum.IntEnum):
        """Integer enum representing the fate of a segment."""

        UNSET = 0  #: Unset.
        CONTINUES = 1  #: Trajectory continues.
        MERGED = 2  #: Trajectory merged (pruned).
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

    @property
    def status_text(self):
        return self.status_names[self.status]

    @property
    def endpoint_type_text(self):
        return self.endpoint_type_names[self.endpoint_type]

    # convenience functions for binning
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
        if weight is not None:
            weight = float(weight)
            if not (0 <= weight <= 1):
                raise ValueError("'weight' must be between 0 and 1")

        self.n_iter = n_iter
        self.seg_id = seg_id
        self._weight = weight  # read-only
        self.endpoint_type = endpoint_type
        self.parent_id = parent_id
        self.wtg_parent_ids = wtg_parent_ids
        self.pcoord = pcoord
        self.status = status
        self.walltime = walltime
        self.cputime = cputime
        self.data = data
        self.initial_state = initial_state
        self.final_state = final_state

    @property
    def n_iter(self):
        return self._n_iter

    @n_iter.setter
    def n_iter(self, value):
        self._n_iter = int(value) if value is not None else None

    @property
    def seg_id(self):
        return self._seg_id

    @seg_id.setter
    def seg_id(self, value):
        self._seg_id = int(value) if value is not None else None

    @property
    def weight(self):
        return self._weight

    @property
    def initpoint_type(self):
        if self.parent_id is None:
            return Segment.InitPoint.UNSET
        elif self.parent_id < 0:
            return Segment.InitPoint.NEWTRAJ
        else:
            return Segment.InitPoint.CONTINUES

    @property
    def endpoint_type(self):
        return self._endpoint_type

    @endpoint_type.setter
    def endpoint_type(self, value):
        self._endpoint_type = Segment.EndPoint(value) if value else Segment.EndPoint.UNSET

    @property
    def parent_id(self):
        return self._parent_id

    @parent_id.setter
    def parent_id(self, value):
        self._parent_id = int(value) if value is not None else None

    @property
    def wtg_parent_ids(self):
        return self._wtg_parent_ids

    @wtg_parent_ids.setter
    def wtg_parent_ids(self, value):
        self._wtg_parent_ids = _IntegerSet(value)

    @property
    def pcoord(self):
        return self._pcoord

    @pcoord.setter
    def pcoord(self, value):
        if value is not None:
            value = np.asarray(value)
            if not np.issubdtype(value.dtype, np.number):
                raise TypeError("scalar type of 'pcoord' must be numeric")
            if value.ndim != 2:
                raise ValueError("'pcoord' must be 2-D array")
            if len(value) < 2:
                raise ValueError("'pcoord' time series length must be at least 2")
        self._pcoord = value

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, value):  # NOTE(Jeff): Default changed from None to UNSET.
        self._status = Segment.Status(value) if value else Segment.Status.UNSET

    @property
    def walltime(self):
        return self._walltime

    @walltime.setter
    def walltime(self, value):
        if value is None:
            value = 0.0
        else:
            value = float(value)
            if value < 0:
                raise ValueError("'walltime' must be nonnegative")
        self._walltime = value

    @property
    def cputime(self):
        return self._cputime

    @cputime.setter
    def cputime(self, value):
        if value is None:
            value = 0.0
        else:
            value = float(value)
            if value < 0:
                raise ValueError("'cputime' must be nonnegative")
        self._cputime = value

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = _AuxiliaryData(value or {})

    @property
    def initial_state(self):
        return self._initial_state

    @initial_state.setter
    def initial_state(self, value):
        if not (isinstance(value, State) or value is None):
            raise TypeError("'initial_state' must be a State object or None")
        self._initial_state = value

    @property
    def final_state(self):
        return self._final_state

    @final_state.setter
    def final_state(self, value):
        if not (isinstance(value, State) or value is None):
            raise TypeError("'final_state' must be a State object or None")
        self._final_state = value

    def __repr__(self):
        return '<%s n_iter=%r, seg_id=%r, weight=%r, parent_id=%r at %s>' % (
            type(self).__name__,
            self.n_iter,
            self.seg_id,
            self.weight,
            self.parent_id,
            hex(id(self)),
        )

    @property
    def initial_state_id(self):  # deprecated: always returns -1
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
        """Return a copy of the segment with modified attributes.

        Parameters
        ----------
        **changes
            Name-value pairs specifying the attributes to modify.

        Returns
        -------
        new_segment : Segment
            Copy of the segment with attributes modified according to `changes`.

        """
        return self.__replace__(**changes)
