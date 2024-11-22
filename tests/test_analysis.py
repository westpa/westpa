import itertools
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from westpa.analysis import Run
from westpa.analysis.core import Iteration, Walker
from westpa.core.h5io import WESTPAH5File
from westpa.core.states import InitialState


@pytest.fixture(scope='module')
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


@pytest.fixture(scope='module')
def run(h5filename) -> Run:
    run = Run.open(h5filename)
    yield run
    run.close()


def test_num_iterations(run):
    assert run.num_iterations == 50
    assert len(run) == run.num_iterations  # __len__ is implemented


def test_iteration(run):
    for n in range(1, run.num_iterations + 1):
        iteration = run.iteration(n)
        assert isinstance(iteration, Iteration)
        assert iteration in run
        assert iteration.run is run
        assert iteration.number == n


def test_iterations(run):
    assert isinstance(run.iterations, list)
    assert run.iterations == list(run)  # __iter__ is implemented
    for iteration, n in zip(run.iterations, range(1, run.num_iterations + 1)):
        assert isinstance(iteration, Iteration)
        assert iteration in run
        assert iteration.run is run
        assert iteration.number == n


def test_num_walkers(run):
    assert run.num_walkers == 9985
    assert sum(iteration.num_walkers for iteration in run) == run.num_walkers
    assert run.num_walkers == run.num_segments  # alias
    for iteration in run:
        assert iteration.num_walkers == iteration.num_segments  # alias


def test_walker(run):
    for iteration in run:
        for i in range(iteration.num_walkers):
            walker = iteration.walker(i)
            assert isinstance(walker, Walker)
            assert walker in run
            assert walker in iteration
            assert walker.run is run
            assert walker.iteration is iteration
            assert walker.index == i


def test_walkers(run):
    for walker in run.walkers:
        assert isinstance(walker, Walker)
        assert walker in run
        assert walker.run is run
    assert len(list(run.walkers)) == run.num_walkers
    for iteration in run:
        for walker, i in zip(iteration.walkers, range(iteration.num_walkers)):
            assert isinstance(walker, Walker)
            assert walker in run
            assert walker in iteration
            assert walker.run is run
            assert walker.iteration is iteration
            assert walker.index == i
    assert sum(len(list(iteration.walkers)) for iteration in run) == run.num_walkers


def test_recycled_walkers(run):
    for walker in run.recycled_walkers:
        assert walker.recycled
        assert walker.pcoords[-1] in walker.iteration.sink
    assert len(list(run.recycled_walkers)) == 39
    for walker1, walker2 in zip(
        run.recycled_walkers,
        itertools.chain(*(iteration.recycled_walkers for iteration in run)),
    ):
        assert walker1 == walker2


def test_initial_walkers(run):
    for walker in run.initial_walkers:
        assert walker.initial
        assert isinstance(walker.parent, InitialState)
    assert len(list(run.initial_walkers)) == 5
    for walker1, walker2 in zip(
        run.initial_walkers,
        itertools.chain(*(iteration.initial_walkers for iteration in run)),
    ):
        assert walker1 == walker2


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
    assert np.all(summary.index == np.arange(1, run.num_iterations + 1))
    assert list(summary.axes[1]) == fields
    for iteration in run:
        summary = iteration.summary
        assert isinstance(summary, pd.Series)
        assert list(summary.index) == fields
        assert summary.name == iteration.number


def test_h5filename(run, h5filename):
    assert run.h5filename == h5filename


def test_h5file(run):
    assert isinstance(run.h5file, WESTPAH5File)
    assert run.h5file.filename == run.h5filename


def test_h5group(run):
    for iteration in run:
        assert isinstance(iteration.h5group, h5py.Group)
        assert iteration.h5group == run.h5file.get_iter_group(iteration.number)
