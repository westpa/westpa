from collections.abc import Sequence
from typing import Protocol

from westpa.core.segment import Segment
from westpa.core.binning import Bin


class Propagator(Protocol):
    """Callback protocol for propagators."""

    def __call__(self, segments: Sequence[Segment]) -> Sequence[Segment]:
        """A propagator is a callable that accepts a sequence of segments,
        "propagates" each segment by reading its ``initial_state`` and
        setting its ``final_state``, and returns the modified segments.
        Besides setting ``final_state``, a propagator may also modify the
        ``pcoord``, ``data``, ``walltime``, or ``cputime`` attributes of a segment.

        To illustrate, the function below is a propagator that keeps the state
        of each walker constant::

            def static_propagator(segments):
                for segment in segments:
                    segment.final_state = segment.initial_state
                return segments

        Optionally, propagators may define a ``block_size`` attribute.
        If defined, ``block_size`` (a positive integer or None) is used by the
        :class:`Simulation` class to determine the maximum number of segments to
        pass to the propagator in a given call (i.e., work manager task).

        See Also
        --------
        SerialPropagator, VectorizedPropagator
            Base classes providing common functionality for implementing propagators.
        AmberPropagator, GROMACSPropagator, OpenMMPropagator
            Molecular dynamics propagators.

        """
        ...


class BinMapper(Protocol):
    """Callback protocol for bin mappers."""

    def __call__(self, segments: Sequence[Segment]) -> Sequence[Bin]:
        """A bin mapper is a callable that accepts a sequence of segments and
        returns a sequence of bins, together containing the input segments.

        To illustrate, the function below is a bin mapper that assigns all the
        segments to a single bin::

            def nop_bin_mapper(segments):
                return [westpa.Bin(segments)]

        """
        ...


class Resampler(Protocol):
    """Callback protocol for resamplers."""

    def __call__(self, bin: Bin, target_count: int) -> Bin:
        """A resampler is a callable that, given a bin and target number of
        walkers for the bin, resamples the walkers in the bin and returns the modified bin.

        """
        ...
