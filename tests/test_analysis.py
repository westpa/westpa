import itertools
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from westpa.analysis import Run
from westpa.analysis.core import Iteration, Walker
from westpa.core.binning import RectilinearBinMapper
from westpa.core.h5io import WESTPAH5File
from westpa.core.states import BasisState, InitialState, TargetState


class TestAnalysis:
    @pytest.fixture
    def h5filename(self) -> str:
        p = Path(__file__).parents[0] / 'refs' / 'west_ref.h5'
        return str(p)

    @pytest.fixture
    def run(self, h5filename) -> Run:
        run = Run.open(h5filename)
        yield run
        run.close()

    def test_open_close(self, h5filename):
        run_obj = Run.open(h5filename)
        assert not run_obj.closed
        run_obj.close()
        assert run_obj.closed

    def test_context_manager(self, h5filename):
        with Run.open(h5filename) as run_obj:
            assert not run_obj.closed
        assert run_obj.closed

    def test_num_iterations(self, run):
        assert run.num_iterations == 50
        assert len(run) == run.num_iterations  # Run.__len__

    def test_iteration(self, run):
        for n in range(1, run.num_iterations + 1):
            iteration = run.iteration(n)
            assert isinstance(iteration, Iteration)
            assert iteration in run
            assert iteration.run is run
            assert iteration.number == n

    def test_iterations(self, run):
        assert isinstance(run.iterations, list)
        assert run.iterations == list(run)  # Run.__iter__, Iteration.__eq__
        for iteration, n in zip(run.iterations, range(1, run.num_iterations + 1)):
            assert isinstance(iteration, Iteration)
            assert iteration in run
            assert iteration.run is run
            assert iteration.number == n

    def test_num_walkers(self, run):
        assert run.num_walkers == 9985
        assert sum(iteration.num_walkers for iteration in run) == run.num_walkers
        assert run.num_walkers == run.num_segments  # alias
        for iteration in run:
            assert iteration.num_walkers == iteration.num_segments  # alias

    def test_walker(self, run):
        iteration = run.iteration(1)
        for i in range(iteration.num_walkers):
            walker = iteration.walker(i)
            assert isinstance(walker, Walker)
            assert walker in run
            assert walker in iteration
            assert walker.run is run
            assert walker.iteration is iteration
            assert walker.index == i

    def test_walkers(self, run):
        for walker in run.walkers:
            assert isinstance(walker, Walker)
            assert walker in run
            assert walker.run is run

        assert len(list(run.walkers)) == run.num_walkers

        for iteration in run:
            assert list(iteration.walkers) == list(iteration)  # Iteration.__iter__, Walker.__eq__
            for i, walker in enumerate(iteration.walkers):
                assert isinstance(walker, Walker)
                assert walker in run
                assert walker in iteration
                assert walker.run is run
                assert walker.iteration is iteration
                assert walker.index == i
        assert sum(len(list(iteration.walkers)) for iteration in run) == run.num_walkers

    def test_parent(self, run):
        for walker in run.iteration(1):
            assert isinstance(walker.parent, InitialState)
        for walker in run.iteration(2):
            assert isinstance(walker.parent, Walker)
            assert walker.parent.iteration.number == 1

    def test_children(self, run):
        for walker in run.iteration(1):
            for child in walker.children:
                assert isinstance(child, Walker)
                assert child.iteration.number == 2

    def test_recycled(self, run):
        for walker in run.iteration(1):
            assert not walker.recycled

    def test_initial(self, run):
        for walker in run.iteration(1):
            assert walker.initial

    def test_recycled_walkers(self, run):
        for walker in run.recycled_walkers:
            assert walker.recycled
        assert len(list(run.recycled_walkers)) == 39
        for walker1, walker2 in zip(
            run.recycled_walkers,
            itertools.chain(*(iteration.recycled_walkers for iteration in run)),
        ):
            assert walker1 == walker2

    def test_initial_walkers(self, run):
        for walker in run.initial_walkers:
            assert walker.initial
            assert isinstance(walker.parent, InitialState)
        assert len(list(run.initial_walkers)) == 5
        for walker1, walker2 in zip(
            run.initial_walkers,
            itertools.chain(*(iteration.initial_walkers for iteration in run)),
        ):
            assert walker1 == walker2

    def test_auxiliary_data(self, run):
        for iteration in run:
            assert iteration.auxiliary_data is None
        for walker in run.iteration(1):
            assert walker.auxiliary_data == {}

    def test_basis_state_summaries(self, run):
        for iteration in run:
            summaries = iteration.basis_state_summaries
            assert isinstance(summaries, pd.DataFrame)
            assert list(summaries.axes[1]) == ['label', 'probability', 'auxref']

    def test_basis_state_pcoords(self, run):
        for iteration in run:
            assert iteration.basis_state_pcoords.ndim == 2

    def test_basis_states(self, run):
        for iteration in run:
            basis_states = iteration.basis_states
            assert isinstance(basis_states, list)
            assert all(isinstance(state, BasisState) for state in basis_states)

    def test_has_target_states(self, run):
        for iteration in run:
            assert iteration.has_target_states

    def test_target_state_summaries(self, run):
        for iteration in run:
            summaries = iteration.target_state_summaries
            assert isinstance(summaries, pd.DataFrame)
            assert list(summaries.axes[1]) == ['label']

    def test_target_state_pcoords(self, run):
        for iteration in run:
            assert iteration.target_state_pcoords.ndim == 2

    def test_target_states(self, run):
        for iteration in run:
            target_states = iteration.target_states
            assert isinstance(target_states, list)
            assert all(isinstance(state, TargetState) for state in target_states)

    def test_sink(self, run):
        for walker in run.recycled_walkers:
            assert walker.pcoords[0] not in walker.iteration.sink
            assert walker.pcoords[-1] in walker.iteration.sink

    def test_summary(self, run):
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
        assert all(summary.index == range(1, run.num_iterations + 1))
        assert list(summary.axes[1]) == fields
        for iteration in run:
            summary = iteration.summary
            assert isinstance(summary, pd.Series)
            assert list(summary.index) == fields
            assert summary.name == iteration.number

    def test_segment_summaries(self, run):
        fields = [
            'weight',
            'parent_id',
            'wtg_n_parents',
            'wtg_offset',
            'cputime',
            'walltime',
            'endpoint_type',
            'status',
        ]
        for iteration in run:
            summaries = iteration.segment_summaries
            assert isinstance(summaries, pd.DataFrame)
            assert len(summaries) == iteration.num_walkers
            assert all(summaries.index == range(iteration.num_walkers))
            assert list(summaries.axes[1]) == fields
        for walker in run.iteration(1):
            summary = walker.segment_summary
            assert isinstance(summary, pd.Series)
            assert list(summary.index) == fields

    def test_h5filename(self, run, h5filename):
        assert run.h5filename == h5filename

    def test_h5file(self, run):
        assert isinstance(run.h5file, WESTPAH5File)
        assert run.h5file.filename == run.h5filename

    def test_h5group(self, run):
        for iteration in run:
            assert isinstance(iteration.h5group, h5py.Group)
            assert iteration.h5group == run.h5file.get_iter_group(iteration.number)

    def test_prev_next(self, run):
        for iteration1, iteration2 in zip(run.iterations[:-1], run.iterations[1:]):
            assert iteration1.next == iteration2
            assert iteration2.prev == iteration1

    def test_pcoords(self, run):
        for iteration in run:
            assert iteration.pcoords.ndim == 3
            assert len(iteration.pcoords) == iteration.num_walkers
        pcoords = run.iteration(1).pcoords
        for walker in run.iteration(1):
            assert np.allclose(walker.pcoords, pcoords[walker.index])

    def test_weights(self, run):
        for iteration in run:
            assert iteration.weights.ndim == 1
            assert len(iteration.weights) == iteration.num_walkers
        weights = run.iteration(1).weights
        for walker in run.iteration(1):
            assert np.isclose(walker.weight, weights[walker.index])

    def test_binning(self, run):
        for iteration in run:
            bin_mapper = iteration.bin_mapper
            if bin_mapper is None:
                assert iteration.bin_target_counts is None
                assert iteration.num_bins == 0
            else:
                assert isinstance(bin_mapper, RectilinearBinMapper)
                assert len(iteration.bin_target_counts) == bin_mapper.nbins
                assert iteration.num_bins == bin_mapper.nbins

    def test_trace(self, run):
        walker = run.iterations[-1].walker(0)
        trace = walker.trace()
        assert len(trace) == walker.iteration.number
        assert all(isinstance(walker, Walker) for walker in trace)
        assert isinstance(trace.initial_state, InitialState)
