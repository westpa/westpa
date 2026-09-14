import logging
import os
import secrets
import time
from abc import ABC, abstractmethod

import numpy as np

logger = logging.getLogger(__name__)


class _PropagatorBase(ABC):
    DEFAULT_SEGMENT_DIR_TEMPLATE = 'traj_segs/{n_iter:06d}/{seg_id:06d}'

    def __init__(
        self,
        root_seed=None,
        bit_generator_type=None,
        segment_dir_template=None,
        block_size=None,
    ):
        self.root_seed = root_seed if root_seed is not None else secrets.randbits(128)
        self.bit_generator_type = bit_generator_type or np.random.PCG64
        self.segment_dir_template = segment_dir_template or self.DEFAULT_SEGMENT_DIR_TEMPLATE
        self.block_size = block_size

    @property
    def root_seed(self):
        """Root seed integer."""
        return self._root_seed

    @root_seed.setter
    def root_seed(self, value):
        if not isinstance(value, int):
            raise TypeError("'root_seed' must be an integer")
        if value < 0:
            raise ValueError("'root_seed' must be non-negative")
        self._root_seed = value
        logger.info(f'root_seed={value}')

    @property
    def bit_generator_type(self):
        """Algorithm for generating random bits."""
        return self._bit_generator_type

    @bit_generator_type.setter
    def bit_generator_type(self, value):
        if not issubclass(value, np.random.BitGenerator):
            raise TypeError("'bit_generator_type' must be a subclass of numpy.random.BitGenerator")
        self._bit_generator_type = value

    @property
    def segment_dir_template(self):
        """Segment output directory template."""
        return self._segment_dir_template

    @segment_dir_template.setter
    def segment_dir_template(self, value):
        self._segment_dir_template = os.path.abspath(value)

    @property
    def block_size(self):
        """Batch size for propagation tasks."""
        return self._block_size

    @block_size.setter
    def block_size(self, value):
        if value is not None:
            if not isinstance(value, int):
                raise TypeError("'block_size' must be an integer or None")
            if value < 1:
                raise ValueError("'block_size' must be positive")
        self._block_size = value

    def _get_rng(self, segment):
        seed = [segment.seg_id, segment.n_iter, self.root_seed]
        bit_generator = self._bit_generator_type(seed)
        return np.random.default_rng(bit_generator)

    def make_segment_dir(self, segment):
        """Create a directory to store output for a given segment, and return
        its path.

        Parameters
        ----------
        segment : Segment
            Segment to create the output directory for.

        Returns
        -------
        segment_dir : str
            Absolute path of the new directory.

        """
        segment_dir = self.segment_dir_template.format(n_iter=segment.n_iter, seg_id=segment.seg_id)
        os.makedirs(segment_dir)
        return segment_dir

    @abstractmethod
    def __call__(self, segments):  # noqa
        ...


class SerialPropagator(_PropagatorBase):
    """Base class for propagators that propagate one segment at a time.
    Subclasses must implement the :meth:`propagate` method.

    Parameters
    ----------
    root_seed : int, optional
        Root seed integer. Must be non-negative. Used to reproducibly seed
        the pseudorandom number generator (PRNG) passed to the
        :meth:`propagate` method. Defaults to ``secrets.randbits(128)``.
    bit_generator_type : type, optional
        Subclass of ``numpy.random.BitGenerator`` specifying the PRNG
        algorithm to use. Defaults to ``numpy.random.PCG64``.
    segment_dir_template : str, optional
        Path template used by the :meth:`make_segment_dir` method. Must
        contain ``{n_iter}`` and ``{seg_id}`` replacement fields. If a
        relative path is provided, it is assumed to be relative to the
        current working directory.
        Defaults to ``'traj_segs/{n_iter:06d}/{seg_id:06d}'``.
    block_size : int, default 1
        Maximum number of segments to pass to the propagator in a given call.

    Attributes
    ----------
    root_seed : int
    bit_generator_type : type
    segment_dir_template : str
    block_size : int or None

    Methods
    -------
    propagate
    make_segment_dir

    Examples
    --------
    Create a propagator that simulates a random walk on the integers:

    >>> import westpa
    >>> class RandomWalkPropagator(westpa.SerialPropagator):
    ...     def __init__(self, p=0.5, steps=1, **kwargs):
    ...         super().__init__(**kwargs)
    ...         self.p = [p, 1 - p]
    ...         self.steps = steps
    ...     def propagate(self, segment, rng):
    ...         delta = rng.choice([-1, 1], p=self.p, size=self.steps).sum()
    ...         final_coord = segment.initial_state.coord + delta
    ...         segment.final_state = westpa.State(coord=final_coord)
    ...         return segment
    ...
    >>> propagator = RandomWalkPropagator(steps=100, root_seed=12345)

    Propagate a test segment:

    >>> segment = westpa.Segment(1, 0)
    >>> segment.initial_state = westpa.State(coord=[0]))
    >>> segment, = propagator([segment])
    >>> segment.final_state.coord
    array([-8])
    >>> segment.walltime  # doctest: +SKIP
    0.00029450003057718277

    """

    def __init__(
        self,
        root_seed=None,
        bit_generator_type=None,
        segment_dir_template=None,
        block_size=1,
    ):
        super().__init__(
            root_seed=root_seed,
            bit_generator_type=bit_generator_type,
            segment_dir_template=segment_dir_template,
            block_size=block_size,
        )

    @abstractmethod
    def propagate(self, segment, rng):
        """Propagate a single trajectory segment.

        Parameters
        ----------
        segment : Segment
            Segment to propagate.
        rng : numpy.random.Generator
            PRNG initialized as follows::

                seed = [segment.seg_id, segment.n_iter, self.root_seed]
                bit_generator = self.bit_generator_type(seed)
                rng = numpy.random.default_rng(bit_generator)

            The choice of seed is based on the
            `sequence of integer seeds <https://numpy.org/doc/stable/reference/random/parallel.html#sequence-of-integer-seeds>`_
            scheme for parallel random number generation.

        Returns
        -------
        segment : Segment
            Propagated segment.

        """
        ...

    def __call__(self, segments):
        for segment in segments:
            start_walltime = time.perf_counter()
            segment = self.propagate(segment, self._get_rng(segment))
            segment.walltime = time.perf_counter() - start_walltime
        return segments


class VectorizedPropagator(_PropagatorBase):
    """Base class for propagators that propagate multiple segments at the same time.
    Subclasses must implement the :meth:`propagate` method.

    Parameters
    ----------
    root_seed : int, optional
        Root seed integer. Must be non-negative. Used to reproducibly seed
        the pseudorandom number generator (PRNG) passed to the
        :meth:`propagate` method. Defaults to ``secrets.randbits(128)``.
    bit_generator_type : type, optional
        Subclass of ``numpy.random.BitGenerator`` specifying the PRNG
        algorithm to use. Defaults to ``numpy.random.PCG64``.
    segment_dir_template : str, optional
        Path template used by the :meth:`make_segment_dir` method. Must
        contain ``{n_iter}`` and ``{seg_id}`` replacement fields. If a
        relative path is provided, it is assumed to be relative to the
        current working directory.
        Defaults to ``'traj_segs/{n_iter:06d}/{seg_id:06d}'``.
    block_size : int, optional
        Maximum number of segments to pass to the propagator in a given call.

    Attributes
    ----------
    root_seed : int
    bit_generator_type : type
    segment_dir_template : str
    block_size : int or None

    Methods
    -------
    propagate
    make_segment_dir

    Examples
    --------
    Create a propagator that simulates a random walk on the integers:

    >>> import westpa
    >>> class RandomWalkPropagator(westpa.VectorizedPropagator):
    ...     def __init__(self, p=0.5, steps=1, **kwargs):
    ...         super().__init__(**kwargs)
    ...         self.p = [p, 1 - p]
    ...         self.steps = steps
    ...     def propagate(self, segments, rng):
    ...         deltas = rng.choice(
    ...             [1, -1], p=self.p, size=(len(segments), self.steps)
    ...         ).sum(axis=1)
    ...         for segment, delta in zip(segments, deltas):
    ...             final_coord = segment.initial_state.coord + delta
    ...             segment.final_state = westpa.State(coord=final_coord)
    ...         return segments
    ...
    >>> propagator = RandomWalkPropagator(steps=100, root_seed=12345)

    Propagate a test segment:

    >>> segment = westpa.Segment(1, 0)
    >>> segment.initial_state = westpa.State(coord=[0]))
    >>> segment, = propagator([segment])
    >>> segment.final_state.coord
    array([-8])

    """

    @abstractmethod
    def propagate(self, segments, rng):
        """Propagate a batch of trajectory segments.

        Parameters
        ----------
        segments : sequence of Segment
            Segments to propagate.
        rng : numpy.random.Generator
            PRNG initialized as follows::

                seed = [segments[0].seg_id, segments[0].n_iter, self.root_seed]
                bit_generator = self.bit_generator_type(seed)
                rng = numpy.random.default_rng(bit_generator)

            The choice of seed is based on the
            `sequence of integer seeds <https://numpy.org/doc/stable/reference/random/parallel.html#sequence-of-integer-seeds>`_
            scheme for parallel random number generation.

        Returns
        -------
        segments : sequence of Segment
            Propagated segments.

        """
        ...

    def __call__(self, segments):
        return self.propagate(segments, rng=self._get_rng(segments[0]))
