from pathlib import Path

import pandas as pd
import pytest

from westpa.analysis import Run


@pytest.fixture
def h5filename() -> str:
    p = Path(__file__).parents[1] / 'refs' / 'west_ref.h5'
    return str(p)


def test_open_close(h5filename):
    run = Run.open(h5filename)
    assert not run.closed
    run.close()
    assert run.closed


def test_context_manager(h5filename):
    with Run.open(h5filename) as run:
        assert not run.closed
    assert run.closed


@pytest.fixture
def run(h5filename) -> Run:
    run = Run.open(h5filename)
    yield run
    run.close()


def test_summary(run):
    df = run.summary
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 50


def test_num_iterations(run):
    assert run.num_iterations == 50


def test_num_walkers(run):
    assert run.num_walkers == 9985


def test_recycled_walkers(run):
    for walker in run.recycled_walkers:
        assert walker.pcoords[-1] in walker.iteration.sink
    assert len(list(run.recycled_walkers)) == 39
