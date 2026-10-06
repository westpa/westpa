from collections.abc import Iterable, Mapping, Sequence
from typing import runtime_checkable, Protocol

from numpy.typing import ArrayLike

from westpa.core.binning import Bin
from westpa.core.segment import Segment
from westpa.core.state import State


@runtime_checkable
class Propagator(Protocol):
    """Protocol for propagators.

    Methods
    -------
    __call__

    See Also
    --------
    SerialPropagator, VectorizedPropagator
        Base classes for propagators.
    AmberPropagator, GROMACSPropagator, OpenMMPropagator
        Molecular dynamics propagators.

    Notes
    -----
    A propagator may optionally define a ``block_size`` attribute. Its value
    (a positive integer or None) indicates to the calling
    :class:`~westpa.Simulation` the maximum number of segments to pass to the
    propagator in a single call. If ``block_size`` is None or undefined,
    no chunking of input is done (i.e., prepared segments are all
    passed in a single call).

    Examples
    --------
    No-op implementation (keeps the state of each walker constant):

    >>> def nop_propagator(segments):
    ...     for segment in segments:
    ...         segment.final_state = segment.initial_state
    ...     return segments
    ...

    """

    def __call__(self, segments: Sequence[Segment]) -> Sequence[Segment]:
        """Propagate a set of segments. Propagating a segment consists of
        reading its ``initial_state`` and setting its ``final_state``.
        Besides setting ``final_state``, a propagator may also modify the
        ``pcoord``, ``data``, ``walltime``, or ``cputime`` attributes.

        Parameters
        ----------
        segments : sequence of Segment
            Segments to propagate.

        Returns
        -------
        segments : sequence of Segment
            Propagated segments.

        """
        ...


@runtime_checkable
class PCoordCalculator(Protocol):
    """Protocol for progress coordinate calculators.

    Methods
    -------
    __call__

    """

    def __call__(self, obj: State | Segment) -> ArrayLike | tuple[ArrayLike, Mapping[str, ArrayLike]]:
        """Return the progress coordinates for a given state or segment.

        Parameters
        ----------
        obj : State or Segment
            State or segment for which to compute progress coordinates.

        Returns
        -------
        pcoord : array_like
            If `obj` is a state, a coordinate point (1-D array). If `obj`
            is a segment, a time series (2-D array) with at least two points,
            corresponding to the initial and final states of the segment.
        data : Mapping[str, array_like], optional
            Dictionary of named arrays to store as auxiliary data.

        Examples
        --------
        Use the raw coordinates as progress coordinates (default behavior):

        >>> import westpa
        >>> def default_pcoord_calculator(obj):
        ...     if isinstance(obj, westpa.State):
        ...         return obj.coord
        ...     else:
        ...         return obj.initial_state.coord, obj.final_state.coord
        ...


        """
        ...


@runtime_checkable
class BinMapper(Protocol):
    """Protocol for bin mappers.

    Attributes
    ----------
    nbins : int
        Number of bins mapped to.
    labels : iterable of str
        Bin labels.

    Methods
    -------
    __call__

    Examples
    --------
    No-op implementation (assigns all segments to a single bin):

    >>> def nop_bin_mapper(segments, coord_index=-1):
    ...     return [0] * len(segments)
    ...
    >>> nop_bin_mapper.nbins = 1
    >>> nop_bin_mapper.labels = ['nop']

    See Also
    --------
    BinMapperBase
        Base class for bin mappers.

    """

    nbins: int
    labels: Iterable[str]

    def __call__(self, segments: Sequence[Segment], coord_index: int = -1) -> ArrayLike:
        """Assign segments to bins.

        Parameters
        ----------
        segments : sequence of Segment
            Segments to be binned.
        coord_index : int, default -1
            Index of the progress coordinate or auxiliary data point to use
            for binning. Defaults to the final point. Bin mappers that do not
            depend on coordinates may ignore this parameter (see, for instance,
            the no-op example above).

        Returns
        -------
        assignments : 1-D array_like of integer type
            Bin assignment for each segment. Array elements must be in
            ``range(self.nbins)``.

        """
        ...


@runtime_checkable
class Resampler(Protocol):
    """Protocol for resamplers.

    Methods
    -------
    __call__

    See Also
    --------
    ResamplerBase
        Base class for resamplers.
    HuberKimResampler
        Default resampler.

    """

    def __call__(self, bin_: Bin, target_count: int) -> Bin:
        """Resample the trajectories in a given bin.

        Parameters
        ----------
        bin_ : Bin
            Bin to be resampled.
        target_count : int
            Target number of trajectories for the bin.

        Returns
        -------
        bin_ : Bin
            Resampled bin.

        """
        ...
