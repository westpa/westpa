import logging
import secrets
from collections.abc import Container, Sequence

import numpy as np

from .state import State

logger = logging.getLogger(__name__)


class Source(Sequence):
    """Set of states from which trajectories are reinitiated upon reaching a sink.

    Parameters
    ----------
    states : State or iterable of State, optional
        One or more source states.
    p : 1-D array_like, optional
        Selection probability of each state. Defaults to a uniform distribution.
    rng : numpy.random.Generator, int, or sequence of int, optional
        Psuedorandom number generator (PRNG) to use for sampling states from
        the source, or a seed for initializing the PRNG. Integer values must
        be nonnegative. Defaults to ``numpy.random.default_rng()``.

    Attributes
    ----------
    p : numpy.ndarray
        State selection probabilities.
    rng : numpy.random.Generator
        PRNG used by the :meth:`random_sample` method.

    Examples
    --------

    >>> import westpa
    >>> states = [westpa.State([0.]), westpa.State([1.])]
    >>> source = westpa.Source(states)
    >>> source
    <Source with 2 states, p=[0.5, 0.5] at 0x30218eea0>

    Nonuniform selection probabilities:

    >>> source = westpa.Source(states, p=[0.7, 0.3])
    >>> source
    <Source with 2 states, p=[0.7, 0.3] at 0x3022720A0>

    Random sampling:

    >>> source.random_sample(3)
    [State(coord=array([0.])),
     State(coord=array([0.])),
     State(coord=array([1.]))]

    """

    def __init__(self, states, p=None, rng=None):
        if isinstance(states, State):
            states = (states,)
        else:
            states = tuple(states)
            if not all(isinstance(item, State) for item in states):
                raise TypeError("'states' must be a State object or an iterable of State objects")

        if p is None:
            p = np.ones(len(states))
        else:
            p = np.asarray(p, dtype=float)
            if len(p) != len(states):
                raise ValueError("length of 'p' must match the number of states")
        p /= p.sum()

        if rng is None:
            seed = secrets.randbits(128)
            rng = np.random.default_rng(seed)
            logger.info(f'rng=default_rng({seed=})')
        else:
            rng = np.random.default_rng(rng)

        self._states = states
        self._p = p
        self._rng = rng

    @property
    def p(self):
        return self._p

    @property
    def rng(self):
        return self._rng

    def random_sample(self, size=None):
        """Return a random sample of states from the source distribution.

        Parameters
        ----------
        size : int, optional
            Sample size. Default is None, in which case a single state
            is returned.

        Returns
        -------
        states : State or list of State
            Sampled states.

        """
        if size is None:
            return self.rng.choice(self._states, p=self.p)
        else:
            return self.rng.choice(self._states, p=self.p, size=size).tolist()

    def __getitem__(self, index):
        return self._states[index]

    def __len__(self):
        return len(self._states)

    def __repr__(self):
        s = '' if len(self) == 1 else 's'
        p = np.array2string(self._p, separator=', ')
        return f'<{type(self).__name__} with {len(self)} state{s}, p={p} at {hex(id(self))}>'


class Sink(Container):
    """Represents a sink (target) region.

    Parameters
    ----------
    indicator : callable
        Indicator function. Must accept a :class:`Segment` object as
        input and return True if the corresponding walker reached
        (hit) the sink, False otherwise.
    label : str, optional
        Descriptive label for the sink.

    Attributes
    ----------
    indicator : callable
        Indicator function.
    label : str
        Sink label.

    Examples
    --------

    Create a sink for walkers with final (-1), first-dimension (0)
    progress coordinate values greater than 1.0:

    >>> import westpa
    >>> sink = westpa.Sink(lambda segment: segment.pcoord[-1, 0] > 1.0)

    Test for membership:

    >>> westpa.Segment(pcoord=[[0.0], [0.9]]) in sink
    False
    >>> westpa.Segment(pcoord=[[0.0], [1.1]]) in sink
    True

    """

    def __init__(self, indicator, label=None):
        self._indicator = indicator
        self._label = str(label) if label else ''

    @property
    def indicator(self):
        return self._indicator

    @property
    def label(self):
        return self._label

    def __contains__(self, segment):
        return self.indicator(segment)

    def __repr__(self):
        args = f'indicator={self.indicator!r}'
        if self.label:
            args += f', label={self.label!r}'
        return type(self).__name__ + '(' + args + ')'
