"""Tests for the new Simulation API (westpa.core.simulation)."""

import pytest
import numpy as np
import westpa

from westpa.core.binning import NopMapper
from westpa.work_managers import SerialWorkManager

# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


class TrivialPropagator(westpa.SerialPropagator):
    """Moves a 1-D coordinate by a fixed delta each step."""

    def __init__(self, delta=0.1, **kwargs):
        super().__init__(**kwargs)
        self.delta = delta

    def propagate(self, segment, rng):
        coord = segment.initial_state.coord
        segment.final_state = westpa.State(coord=coord + self.delta)
        return segment


class FailingPropagator(westpa.SerialPropagator):
    """Always raises an exception."""

    def propagate(self, segment, rng):
        raise RuntimeError("deliberate test failure")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def datafile(tmp_path):
    return str(tmp_path / "west.h5")


@pytest.fixture
def propagator():
    return TrivialPropagator()


@pytest.fixture
def sim(datafile, propagator):
    return westpa.Simulation(
        datafile=datafile,
        propagator=propagator,
    )


# ---------------------------------------------------------------------------
# Constructor tests
# ---------------------------------------------------------------------------


class TestSimulationConstructor:
    def test_basic_construction(self, datafile, propagator):
        sim = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
        )
        assert sim.datafile == datafile
        assert sim.propagator is propagator

    def test_default_bin_mapper_is_nop(self, sim):
        assert isinstance(sim.bin_mapper, NopMapper)

    def test_default_bin_target_counts(self, sim):
        np.testing.assert_array_equal(sim.bin_target_counts, [1])

    def test_default_resampler_is_huber_kim(self, sim):
        assert isinstance(sim.resampler, westpa.HuberKimResampler)

    def test_default_work_manager_is_serial(self, sim):
        assert isinstance(sim.work_manager, SerialWorkManager)

    def test_default_source_is_none(self, sim):
        assert sim.source is None

    def test_default_sinks_is_empty(self, sim):
        assert len(sim.sinks) == 0

    def test_default_istate_generator_is_none(self, sim):
        assert sim.istate_generator is None

    def test_invalid_propagator_type(self, datafile):
        with pytest.raises(TypeError, match="'propagator' must be callable"):
            westpa.Simulation(
                datafile=datafile,
                propagator="not_a_propagator",
            )

    def test_invalid_pcoord_calculator_not_callable(self, datafile, propagator):
        with pytest.raises(TypeError, match="'pcoord_calculator' must be callable or None"):
            westpa.Simulation(
                datafile=datafile,
                propagator=propagator,
                pcoord_calculator="not_callable",
            )

    def test_invalid_resampler_type(self, datafile, propagator):
        with pytest.raises(TypeError, match="'resampler' must be callable"):
            westpa.Simulation(
                datafile=datafile,
                propagator=propagator,
                resampler="not_a_resampler",
            )

    def test_invalid_work_manager_type(self, datafile, propagator):
        with pytest.raises(TypeError, match="'work_manager' must be a WorkManager object"):
            westpa.Simulation(
                datafile=datafile,
                propagator=propagator,
                work_manager="not_a_work_manager",
            )

    def test_source_without_sink_raises(self, datafile, propagator):
        source = westpa.Source(westpa.State(coord=[0.0]))
        with pytest.raises(ValueError, match="'source' and 'sink' must be provided together"):
            westpa.Simulation(
                datafile=datafile,
                propagator=propagator,
                source=source,
            )

    def test_sink_without_source_raises(self, datafile, propagator):
        sink = westpa.Sink(lambda seg: seg.pcoord[-1, 0] > 1.0)
        with pytest.raises(ValueError, match="'source' and 'sink' must be provided together"):
            westpa.Simulation(
                datafile=datafile,
                propagator=propagator,
                sink=sink,
            )

    def test_source_and_sink_together(self, datafile, propagator):
        source = westpa.Source(westpa.State(coord=[0.0]))
        sink = westpa.Sink(lambda seg: seg.pcoord[-1, 0] > 1.0)
        sim = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
            source=source,
            sink=sink,
        )
        assert sim.source is source
        assert sim.sinks[0] is sink

    def test_invalid_istate_generator_type(self, datafile, propagator):
        source = westpa.Source(westpa.State(coord=[0.0]))
        sink = westpa.Sink(lambda seg: seg.pcoord[-1, 0] > 1.0)
        with pytest.raises(TypeError, match="'istate_generator' must be callable"):
            westpa.Simulation(
                datafile=datafile,
                propagator=propagator,
                source=source,
                sink=sink,
                istate_generator="not_callable",
            )

    def test_custom_bin_mapper_and_target_counts(self, datafile, propagator):
        mapper = westpa.RectilinearBinMapper([[-np.inf, 0.5, np.inf]])  # 2 bins
        sim = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
            bin_mapper=mapper,
            bin_target_counts=3,
        )
        assert sim.bin_mapper is mapper
        np.testing.assert_array_equal(sim.bin_target_counts, [3, 3])


# ---------------------------------------------------------------------------
# configure_recycling + disable_recycling tests
# ---------------------------------------------------------------------------


class TestConfigureRecycling:
    def test_valid_config(self, sim):
        source = westpa.Source(westpa.State([0.0]))
        sink = westpa.Sink(lambda seg: seg.pcoord[-1, 0] > 1.0)

        sim.configure_recycling(source, sink)
        assert sim.source is source
        assert sim.sinks[0] is sink

        sim.disable_recycling()
        assert sim.source is None
        assert sim.sinks == ()

    def test_invalid_source_type(self, sim):
        sink = westpa.Sink(lambda seg: seg.pcoord[-1, 0] > 1.0)
        with pytest.raises(TypeError, match="'source' must be a Source object"):
            sim.configure_recycling("not_a_source", sink)

    def test_invalid_sink_type(self, sim):
        source = westpa.Source(westpa.State([0.0]))
        with pytest.raises(TypeError, match="'sink' must be a Sink object or an iterable of Sink objects"):
            sim.configure_recycling(source, "not_a_sink")


# ---------------------------------------------------------------------------
# initialize tests
# ---------------------------------------------------------------------------


class TestInitialize:
    def test_initialize_single_state(self, sim, tmp_path):
        sim.initialize(westpa.State(coord=[0.5]))
        assert (tmp_path / "west.h5").exists()
        assert len(sim.segments) == 1

    def test_initialize_multiple_states(self, sim):
        states = [westpa.State(coord=[float(i)]) for i in range(4)]
        sim.initialize(states)
        assert len(sim.segments) == 4

    def test_initialize_uniform_weights(self, sim):
        states = [westpa.State(coord=[float(i)]) for i in range(4)]
        sim.initialize(states)
        weights = [seg.weight for seg in sim.segments]
        assert all(pytest.approx(w) == 0.25 for w in weights)

    def test_initialize_custom_weights_normalized(self, sim):
        states = [westpa.State(coord=[float(i)]) for i in range(3)]
        sim.initialize(states, weights=[1, 2, 1])
        weights = [seg.weight for seg in sim.segments]
        assert pytest.approx(sum(weights)) == 1.0
        # The segment with weight=2 gets 2/4 = 0.5
        assert pytest.approx(max(weights)) == 0.5

    def test_initialize_runtime_raises(self, sim):
        state = westpa.State(coord=[0.5])
        sim.initialize(state)
        with pytest.raises(RuntimeError, match="can't initialize the simulation"):
            sim.initialize(state)

    def test_initialize_weights_length_mismatch(self, sim):
        states = [westpa.State(coord=[float(i)]) for i in range(3)]
        with pytest.raises(ValueError, match="length of 'weights' must match"):
            sim.initialize(states, weights=[0.5, 0.5])

    def test_initialize_segments_are_prepared(self, sim):
        states = [westpa.State(coord=[float(i)]) for i in range(2)]
        sim.initialize(states)
        for seg in sim.segments:
            assert seg.status == westpa.Segment.Status.PREPARED

    def test_initialize_segment_initial_states_match(self, sim):
        states = [westpa.State(coord=[0.0]), westpa.State(coord=[1.0])]
        sim.initialize(states)
        init_coords = {tuple(seg.initial_state.coord) for seg in sim.segments}
        assert (0.0,) in init_coords
        assert (1.0,) in init_coords

    def test_initialize_sets_iteration_to_one(self, sim):
        sim.initialize(westpa.State(coord=[0.5]))
        assert sim.n_iter == 1


# ---------------------------------------------------------------------------
# run integration tests
# ---------------------------------------------------------------------------


class TestRun:
    def test_run_one_iteration(self, sim):
        sim.initialize(westpa.State(coord=[0.0]))
        sim.run(n_iters=1)
        assert len(sim.segments) == 1  # NopMapper with target_count=1

    def test_run_multiple_iterations(self, sim):
        sim.initialize(westpa.State(coord=[0.0]))
        sim.run(n_iters=3)
        assert sim.n_iter == 4  # started at 1, ran 3 iterations

    def test_run_preserves_total_probability(self, sim):
        states = [westpa.State(coord=[float(i) * 0.1]) for i in range(4)]
        sim.initialize(states)
        sim.run(n_iters=2)
        total_weight = sum(seg.weight for seg in sim.segments)
        assert pytest.approx(total_weight) == 1.0

    def test_run_with_rectilinear_bin_mapper(self, datafile, propagator):
        sim = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
            bin_mapper=westpa.RectilinearBinMapper([[-np.inf, 0.5, np.inf]]),
            bin_target_counts=2,
        )
        states = [westpa.State(coord=[0.1]), westpa.State(coord=[0.2])]
        sim.initialize(states)
        sim.run(n_iters=2)
        assert sim.n_iter == 3

    def test_failing_propagator(self, datafile):
        sim = westpa.Simulation(
            datafile=datafile,
            propagator=FailingPropagator(),
        )
        sim.initialize(westpa.State(coord=[0.0]))
        with pytest.raises(RuntimeError):
            sim.run(n_iters=1)

    def test_run_with_source_and_sink(self, datafile, propagator):
        """Walkers that reach the sink should be recycled to the source."""
        source = westpa.Source(westpa.State(coord=[0.0]))
        # Sink: any segment whose final pcoord > 0.5 is recycled
        sink = westpa.Sink(lambda seg: seg.pcoord[-1, 0] > 0.5)

        sim = westpa.Simulation(
            datafile=datafile,
            propagator=TrivialPropagator(delta=1.0),  # large step → always sinks
            source=source,
            sink=sink,
        )
        states = [westpa.State(coord=[0.0])]
        sim.initialize(states)

        sim.run(n_iters=2)
        total_weight = sum(seg.weight for seg in sim.segments)
        assert pytest.approx(total_weight) == 1.0

        for segment in sim.segments:
            assert segment.initial_state is None
            assert segment.status == segment.Status.UNSET

    def test_continue_run(self, datafile, propagator):
        sim = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
        )
        sim.initialize(westpa.State(coord=[0.0]))
        sim.run(n_iters=2)
        del sim

        # continue run from last checkpoint in 'datafile'
        sim2 = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
        )
        with pytest.raises(RuntimeError, match="already initialized"):
            sim2.initialize(westpa.State(coord=[0.0]))
        sim2.run(n_iters=2)
        assert sim2.n_iter == 5


class TestPCoordCalculator:

    def test_default(self, sim):
        sim.initialize(westpa.State(coord=[0.0]))
        sim.run(n_iters=1)

        with westpa.TrajectoryTree(sim.datafile) as trajtree:
            segment = trajtree.get_segment(1, 0)

        assert segment.pcoord.shape == (2, 1)
        assert np.allclose(segment.pcoord[0], segment.initial_state.coord)
        assert np.allclose(segment.pcoord[-1], segment.final_state.coord)

    def test_without_auxdata(self, datafile, propagator):
        def pcoord_calculator(segment, parent=None):
            initial_pcoord = parent.pcoord[-1] if parent else segment.initial_state.coord
            final_pcoord = segment.final_state.coord
            return np.stack((initial_pcoord, final_pcoord))

        sim = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
            pcoord_calculator=pcoord_calculator,
        )
        sim.initialize(westpa.State(coord=[0.0]))
        sim.run(n_iters=2)

        with westpa.TrajectoryTree(sim.datafile) as trajtree:
            segment = trajtree.get_segment(1, 0)

        assert segment.pcoord.shape == (2, 1)
        assert np.allclose(segment.pcoord[0], segment.initial_state.coord)
        assert np.allclose(segment.pcoord[-1], segment.final_state.coord)

    def test_with_auxdata(self, datafile, propagator):
        def pcoord_calculator(segment, parent=None):
            initial_pcoord = parent.pcoord[-1] if parent else segment.initial_state.coord
            final_pcoord = segment.final_state.coord
            auxdata = {'a': [1, 2, 3]}
            return np.stack((initial_pcoord, final_pcoord)), auxdata

        sim = westpa.Simulation(
            datafile=datafile,
            propagator=propagator,
            pcoord_calculator=pcoord_calculator,
        )
        sim.initialize(westpa.State(coord=[0.0]))
        sim.run(n_iters=2)

        with westpa.TrajectoryTree(sim.datafile, load_auxdata=True) as trajtree:
            segment = trajtree.get_segment(1, 0)

        assert segment.pcoord.shape == (2, 1)
        assert np.allclose(segment.pcoord[0], segment.initial_state.coord)
        assert np.allclose(segment.pcoord[-1], segment.final_state.coord)
        assert segment.data['a'].tolist() == [1, 2, 3]
