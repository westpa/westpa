import itertools

import h5py
import numpy as np
import pandas as pd

from westpa.analysis import Run
from westpa.analysis.core import Iteration, Walker
from westpa.core.binning import RectilinearBinMapper
from westpa.core.h5io import WESTPAH5File
from westpa.core.states import BasisState, InitialState, TargetState


class TestAnalysis:
    def test_open_close(self, ref_analysis):
        run = Run.open(self.h5_filepath)
        assert not run.closed
        run.close()
        assert run.closed

    def test_context_manager(self, ref_analysis):
        with Run.open(self.h5_filepath) as run:
            assert not run.closed
        assert run.closed

    def test_num_iterations(self, ref_analysis):
        assert self.run.num_iterations == 50
        assert len(self.run) == self.run.num_iterations  # Run.__len__

    def test_iteration(self, ref_analysis):
        for n in range(1, self.run.num_iterations + 1):
            iteration = self.run.iteration(n)
            assert isinstance(iteration, Iteration)
            assert iteration in self.run
            assert iteration.run is self.run
            assert iteration.number == n

    def test_iterations(self, ref_analysis):
        assert isinstance(self.run.iterations, list)
        assert self.run.iterations == list(self.run)  # Run.__iter__, Iteration.__eq__
        for iteration, n in zip(self.run.iterations, range(1, self.run.num_iterations + 1)):
            assert isinstance(iteration, Iteration)
            assert iteration in self.run
            assert iteration.run is self.run
            assert iteration.number == n

    def test_num_walkers(self, ref_analysis):
        assert self.run.num_walkers == 9985
        assert sum(iteration.num_walkers for iteration in self.run) == self.run.num_walkers
        assert self.run.num_walkers == self.run.num_segments  # alias
        for iteration in self.run:
            assert iteration.num_walkers == iteration.num_segments  # alias

    def test_walker(self, ref_analysis):
        iteration = self.run.iteration(1)
        for i in range(iteration.num_walkers):
            walker = iteration.walker(i)
            assert isinstance(walker, Walker)
            assert walker in self.run
            assert walker in iteration
            assert walker.run is self.run
            assert walker.iteration is iteration
            assert walker.index == i

    def test_walkers(self, ref_analysis):
        for walker in self.run.walkers:
            assert isinstance(walker, Walker)
            assert walker in self.run
            assert walker.run is self.run

        assert len(list(self.run.walkers)) == self.run.num_walkers

        for iteration in self.run:
            assert list(iteration.walkers) == list(iteration)  # Iteration.__iter__, Walker.__eq__
            for i, walker in enumerate(iteration.walkers):
                assert isinstance(walker, Walker)
                assert walker in self.run
                assert walker in iteration
                assert walker.run is self.run
                assert walker.iteration is iteration
                assert walker.index == i
        assert sum(len(list(iteration.walkers)) for iteration in self.run) == self.run.num_walkers

    def test_parent(self, ref_analysis):
        for walker in self.run.iteration(1):
            assert isinstance(walker.parent, InitialState)
        for walker in self.run.iteration(2):
            assert isinstance(walker.parent, Walker)
            assert walker.parent.iteration.number == 1

    def test_children(self, ref_analysis):
        for walker in self.run.iteration(1):
            for child in walker.children:
                assert isinstance(child, Walker)
                assert child.iteration.number == 2

    def test_recycled(self, ref_analysis):
        for walker in self.run.iteration(1):
            assert not walker.recycled

    def test_initial(self, ref_analysis):
        for walker in self.run.iteration(1):
            assert walker.initial

    def test_recycled_walkers(self, ref_analysis):
        for walker in self.run.recycled_walkers:
            assert walker.recycled
        assert len(list(self.run.recycled_walkers)) == 39
        for walker1, walker2 in zip(
            self.run.recycled_walkers,
            itertools.chain(*(iteration.recycled_walkers for iteration in self.run)),
        ):
            assert walker1 == walker2

    def test_initial_walkers(self, ref_analysis):
        for walker in self.run.initial_walkers:
            assert walker.initial
            assert isinstance(walker.parent, InitialState)
        assert len(list(self.run.initial_walkers)) == 5
        for walker1, walker2 in zip(
            self.run.initial_walkers,
            itertools.chain(*(iteration.initial_walkers for iteration in self.run)),
        ):
            assert walker1 == walker2

    def test_auxiliary_data(self, ref_analysis):
        for iteration in self.run:
            assert iteration.auxiliary_data is None
        for walker in self.run.iteration(1):
            assert walker.auxiliary_data == {}

    def test_basis_state_summaries(self, ref_analysis):
        for iteration in self.run:
            summaries = iteration.basis_state_summaries
            assert isinstance(summaries, pd.DataFrame)
            assert list(summaries.axes[1]) == ['label', 'probability', 'auxref']

    def test_basis_state_pcoords(self, ref_analysis):
        for iteration in self.run:
            assert iteration.basis_state_pcoords.ndim == 2

    def test_basis_states(self, ref_analysis):
        for iteration in self.run:
            basis_states = iteration.basis_states
            assert isinstance(basis_states, list)
            assert all(isinstance(state, BasisState) for state in basis_states)

    def test_has_target_states(self, ref_analysis):
        for iteration in self.run:
            assert iteration.has_target_states

    def test_target_state_summaries(self, ref_analysis):
        for iteration in self.run:
            summaries = iteration.target_state_summaries
            assert isinstance(summaries, pd.DataFrame)
            assert list(summaries.axes[1]) == ['label']

    def test_target_state_pcoords(self, ref_analysis):
        for iteration in self.run:
            assert iteration.target_state_pcoords.ndim == 2

    def test_target_states(self, ref_analysis):
        for iteration in self.run:
            target_states = iteration.target_states
            assert isinstance(target_states, list)
            assert all(isinstance(state, TargetState) for state in target_states)

    def test_sink(self, ref_analysis):
        for walker in self.run.recycled_walkers:
            assert walker.pcoords[0] not in walker.iteration.sink
            assert walker.pcoords[-1] in walker.iteration.sink

    def test_summary(self, ref_analysis):
        fields = [
            'n_particles',
            'min_bin_prob',
            'max_bin_prob',
            'min_seg_prob',
            'max_seg_prob',
            'cputime',
            'walltime',
        ]
        summary = self.run.summary
        assert isinstance(summary, pd.DataFrame)
        assert len(summary) == len(self.run)
        assert all(summary.index == range(1, self.run.num_iterations + 1))
        assert list(summary.axes[1]) == fields
        for iteration in self.run:
            summary = iteration.summary
            assert isinstance(summary, pd.Series)
            assert list(summary.index) == fields
            assert summary.name == iteration.number

    def test_segment_summaries(self, ref_analysis):
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
        for iteration in self.run:
            summaries = iteration.segment_summaries
            assert isinstance(summaries, pd.DataFrame)
            assert len(summaries) == iteration.num_walkers
            assert all(summaries.index == range(iteration.num_walkers))
            assert list(summaries.axes[1]) == fields
        for walker in self.run.iteration(1):
            summary = walker.segment_summary
            assert isinstance(summary, pd.Series)
            assert list(summary.index) == fields

    def test_h5filename(self, ref_analysis):
        assert self.run.h5filename == self.h5_filepath

    def test_h5file(self, ref_analysis):
        assert isinstance(self.run.h5file, WESTPAH5File)
        assert self.run.h5file.filename == self.run.h5filename

    def test_h5group(self, ref_analysis):
        for iteration in self.run:
            assert isinstance(iteration.h5group, h5py.Group)
            assert iteration.h5group == self.run.h5file.get_iter_group(iteration.number)

    def test_prev_next(self, ref_analysis):
        for iteration1, iteration2 in zip(self.run.iterations[:-1], self.run.iterations[1:]):
            assert iteration1.next == iteration2
            assert iteration2.prev == iteration1

    def test_pcoords(self, ref_analysis):
        for iteration in self.run:
            assert iteration.pcoords.ndim == 3
            assert len(iteration.pcoords) == iteration.num_walkers
        pcoords = self.run.iteration(1).pcoords
        for walker in self.run.iteration(1):
            assert np.allclose(walker.pcoords, pcoords[walker.index])

    def test_weights(self, ref_analysis):
        for iteration in self.run:
            assert iteration.weights.ndim == 1
            assert len(iteration.weights) == iteration.num_walkers
        weights = self.run.iteration(1).weights
        for walker in self.run.iteration(1):
            assert np.isclose(walker.weight, weights[walker.index])

    def test_binning(self, ref_analysis):
        for iteration in self.run:
            bin_mapper = iteration.bin_mapper
            if bin_mapper is None:
                assert iteration.bin_target_counts is None
                assert iteration.num_bins == 0
            else:
                assert isinstance(bin_mapper, RectilinearBinMapper)
                assert len(iteration.bin_target_counts) == bin_mapper.nbins
                assert iteration.num_bins == bin_mapper.nbins

    def test_trace(self, ref_analysis):
        walker = self.run.iterations[-1].walker(0)
        trace = walker.trace()
        assert len(trace) == walker.iteration.number
        assert all(isinstance(walker, Walker) for walker in trace)
        assert isinstance(trace.initial_state, InitialState)
