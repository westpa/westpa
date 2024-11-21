from pathlib import Path

import pandas as pd
import pytest

from westpa.analysis import Run
from westpa.analysis.core import Iteration, Walker


@pytest.fixture
def h5filename() -> str:
    p = Path(__file__).parents[0] / 'refs' / 'west_ref.h5'
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


def test_num_iterations(run):
    assert run.num_iterations == 50
    assert len(run) == 50


def test_summary(run):
    fields = [
        'n_particles',
        'min_bin_prob',
        'max_bin_prob',
        'min_seg_prob',
        'max_seg_prob',
        'cputime',
        'walltime',
    ]
    summary = run.summary
    assert isinstance(summary, pd.DataFrame)
    assert len(summary) == len(run)
    assert list(summary.axes[1]) == fields
    for iteration in run:
        summary = iteration.summary
        assert isinstance(summary, pd.Series)
        assert list(summary.index) == fields


def test_num_walkers(run):
    assert run.num_walkers == 9985
    assert sum(iteration.num_walkers for iteration in run) == run.num_walkers


def test_iteration(run):
    for n in range(1, run.num_iterations + 1):
        iteration = run.iteration(n)
        assert isinstance(iteration, Iteration)
        assert iteration.run is run
        assert iteration.number == n


def test_iter(run):
    assert list(run) == [run.iteration(n) for n in range(1, run.num_iterations + 1)]


def test_iterations(run):
    assert run.iterations == list(run)


def test_recycled_walkers(run):
    for walker in run.recycled_walkers:
        assert walker.pcoords[-1] in walker.iteration.sink
    assert len(list(run.recycled_walkers)) == 39
    for iteration in run:
        for walker in iteration.recycled_walkers:
            assert walker.pcoords[-1] in iteration.sink
