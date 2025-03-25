import pytest

from westpa.core._sink import Sink


def test_from_string():
    with pytest.raises(ValueError, match="invalid syntax: expected '<variables> : <predicate>'"):
        Sink.from_string('x')
    with pytest.raises(ValueError, match="array variable 'x' must be subscripted"):
        Sink.from_string('x : x > 0')
    with pytest.raises(ValueError, match="unpacked variable 'x' may not be subscripted"):
        Sink.from_string('(x, y) : x[0] > 0')
    with pytest.raises(ValueError, match='invalid comparison: x == 0'):
        Sink.from_string('(x, y) : x == 0 and y > 0')


def test_array_variable():
    sink = Sink.from_string('x : x[0]**2 + x[1]**2 < 1 and x[0] > 0')
    assert (0.5, 0.5) in sink
    assert (-0.5, 0.5) not in sink
    assert (1.0, 1.0) not in sink
    assert (0.5, 0.5, 0.0, 0.0) in sink  # non-indexed entries are ignored
    assert (1.0, 1.0, 0.0, 0.0) not in sink
    with pytest.raises(TypeError, match='not subscriptable'):
        1.0 in sink
    with pytest.raises(IndexError, match='index out of range'):
        (1.0,) in sink


def test_unpacked_variables():
    sink = Sink.from_string('(x, y) : x**2 + y**2 < 1 and x > 0')
    assert (0.5, 0.5) in sink
    assert (-0.5, 0.5) not in sink
    assert (1.0, 1.0) not in sink
    with pytest.raises(TypeError, match='cannot unpack non-iterable'):
        1.0 in sink
    with pytest.raises(ValueError, match='not enough values to unpack'):
        (1.0,) in sink
    with pytest.raises(ValueError, match='too many values to unpack'):
        (1.0, 1.0, 1.0) in sink
