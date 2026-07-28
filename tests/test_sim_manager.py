import argparse
import glob
import os
from unittest.mock import MagicMock
from shutil import copy
import pytest

import numpy as np

import westpa
from westpa.core.binning.assign import RectilinearBinMapper
from westpa.core._rc import WESTRC
from westpa.core.segment import Segment
from westpa.core.states import BasisState
from westpa.core.sim_manager import PropagationError

REFERENCE_PATH = os.path.join(os.path.dirname(__file__), 'refs', 'odld')


def dummy_callback_one(self):
    pass


def copy_ref(dest_dir):
    for filename in glob.glob(os.path.join(REFERENCE_PATH, '*.*')):
        copy(filename, dest_dir)


@pytest.fixture(autouse=True)
def sim_manager_setup(request, tmp_path):
    def create_segment(seg_id, shape, init_pcoord, final_pcoord, weight=1.0):
        seg = Segment(
            n_iter=1,
            seg_id=1123,
            weight=weight,
            parent_id=1,
            pcoord=np.zeros(shape),
        )

        seg.pcoord[0] = init_pcoord
        seg.pcoord[-1] = final_pcoord

        return seg

    request.cls.rc = rc = WESTRC()
    parser = argparse.ArgumentParser()
    rc.add_args(parser)

    os.chdir(tmp_path)
    copy_ref(tmp_path)  # Copy all odld reference files over
    os.environ['WEST_SIM_ROOT'] = str(tmp_path)

    config_file_name = os.path.join(tmp_path, request.param)
    args = parser.parse_args(['-r={}'.format(config_file_name)])
    rc.process_args(args)

    request.cls.sim_manager = sim_manager = rc.get_sim_manager()

    request.cls.test_dir = tmp_path
    request.cls.hdf5 = os.path.join("west.h5")
    request.cls.basis_states = [BasisState(label="label", probability=1.0)]
    shape = sim_manager.system.new_pcoord_array().shape
    segments = [create_segment(i, shape, 0.0, 1.5, weight=0.125) for i in range(4)] + [
        create_segment(i + 4, shape, 1.5, 0.5, weight=0.125) for i in range(4)
    ]

    sim_manager.we_driver.new_iteration()
    sim_manager.we_driver.assign(segments)
    sim_manager.we_driver.construct_next()
    sim_manager.segments = {segment.seg_id: segment for segment in segments}
    sim_manager.incomplete_segments = sim_manager.segments
    sim_manager.current_iter_istates = sim_manager.segments
    sim_manager.completed_segments = sim_manager.segments
    sim_manager.report_bin_statistics = MagicMock(return_value=True)

    request.cls.work_manager = work_manager = rc.get_work_manager()
    work_manager.running = True

    data = rc.get_data_manager()
    data.we_h5filename = request.cls.hdf5
    data.prepare_backing()
    data.create_ibstate_group([])
    data.create_initial_states(1)
    data.save_target_states([])
    data.update_segments = MagicMock(return_value=None)

    n_iter = 0
    it_name = data.iter_group_name(n_iter)
    for group in ["seg_index", "parents", "ibstates", "pcoord"]:
        data.we_h5file.create_group(it_name + "/" + group)
    data.get_new_weight_data = MagicMock(return_value=None)
    data.get_segments = MagicMock(return_value=segments)
    sim_manager.we_driver.rc.get_data_manager = MagicMock(return_value=data)
    sim_manager.n_iter = n_iter

    yield

    del rc
    del sim_manager
    del data
    del work_manager
    del os.environ['WEST_SIM_ROOT']


@pytest.mark.parametrize(
    "sim_manager_setup",
    ['west.cfg', 'west_mab.cfg', 'west_binless.cfg'],
    indirect=['sim_manager_setup'],
    ids=['default', 'MABSimManager', 'BinlessSimManager'],
)
class TestSimManager:
    def dummy_callback_one(self):
        system = self.sim_manager.system
        bounds = [0.0, 1.0, 2.0, 3.0]
        system.bin_mapper = RectilinearBinMapper([bounds])

    def dummy_callback_two(self):
        system = self.sim_manager.system
        bounds = [0.0, 1.0, 2.0, 5.0]
        system.bin_mapper = RectilinearBinMapper([bounds])

    def test_sim_manager(self):
        assert self.sim_manager.n_propagated == 0
        assert len(self.sim_manager._callback_table) == 0

    def test_register_callback(self):
        hook = self.sim_manager.prepare_new_iteration

        self.sim_manager.register_callback(hook, self.dummy_callback_one, 3)
        self.sim_manager.register_callback(hook, self.dummy_callback_two, 0)
        self.sim_manager.register_callback(
            hook, dummy_callback_one, 3
        )  # Same name and priority, but different function, should be added
        self.sim_manager.register_callback(
            hook, self.dummy_callback_one, 2
        )  # Duplicate should never be added, even with different priority
        assert hook in self.sim_manager._callback_table

        callbacks = self.sim_manager._callback_table.get(hook, [])

        assert len(callbacks) == 3  # Make sure only 3 added.

        assert (3, self.dummy_callback_one.__name__, self.dummy_callback_one) in callbacks  # noqa
        assert (0, self.dummy_callback_two.__name__, self.dummy_callback_two) in callbacks  # noqa
        assert (3, dummy_callback_one.__name__, dummy_callback_one) in callbacks  # noqa

    def test_invoke_callback(self):
        hook = self.sim_manager.prepare_new_iteration

        self.sim_manager.register_callback(hook, self.dummy_callback_one, 3)
        self.sim_manager.register_callback(hook, self.dummy_callback_two, 0)

        self.sim_manager.invoke_callbacks(hook)

        system = self.sim_manager.system
        assert np.all(system.bin_mapper.boundaries == np.array([0.0, 1.0, 2.0, 3.0]))  # noqa

    def test_process_config(self):
        self.sim_manager.process_config()
        assert self.sim_manager.do_gen_istates
        assert self.sim_manager.propagator_block_size == 10000
        assert not self.sim_manager.save_transition_matrices
        assert self.sim_manager.max_run_walltime == 10800
        assert self.sim_manager.max_total_iterations == 100

    def test_load_plugins(self):
        self.sim_manager.load_plugins()

    def test_report_bin_statistics(self):
        self.sim_manager.report_bin_statistics([0.0, 1.0, 2.0, 5.0])

    def test_get_bstate_pcoords(self, monkeypatch):
        with monkeypatch.context() as m:
            m.setattr(westpa, 'rc', self.rc)
            self.sim_manager.get_bstate_pcoords(self.basis_states)

    def test_report_basis_states(self):
        self.sim_manager.report_basis_states(self.basis_states)

    def test_report_target_states(self):
        self.sim_manager.report_target_states(self.basis_states)

    @pytest.mark.skip(reason="Cannot currently test WESimManager.initialize_simulation()")
    def test_initialize_simulation(self):
        # TODO: determine how to test self.simulation_manager.initialize_simulation()
        pass

    def test_prepare_iteration(self, monkeypatch):
        with monkeypatch.context() as m:
            m.setattr(westpa, 'rc', self.rc)
            self.sim_manager.prepare_new_iteration()
            self.sim_manager.prepare_iteration()

    def test_finalize_iteration(self, monkeypatch):
        with monkeypatch.context() as m:
            m.setattr(westpa, 'rc', self.rc)
            self.sim_manager.finalize_iteration()

    def test_get_istate_futures(self):
        self.sim_manager.get_istate_futures()

    def test_propagate(self):
        westpa.core.states.pare_basis_initial_states = MagicMock(return_value=([], []))
        self.sim_manager.propagate

    def test_save_bin_data(self):
        self.sim_manager.save_bin_data()

    def test_check_propagation(self, monkeypatch):
        with monkeypatch.context() as m:
            m.setattr(westpa, 'rc', self.rc)

            with pytest.raises(PropagationError):
                self.sim_manager.check_propagation()

    def test_run_we(self):
        self.sim_manager.run_we()

    def test_run(self, monkeypatch):
        with monkeypatch.context() as m:
            m.setattr(westpa, 'rc', self.rc)
            self.sim_manager.run()

    def test_prepare_run(self):
        self.sim_manager.prepare_run()

    def test_finalize_run(self):
        self.sim_manager.finalize_run()

    def test_pre_propagation(self):
        self.sim_manager.pre_propagation()

    def test_post_propagation(self):
        self.sim_manager.post_propagation()

    def test_pre_we(self):
        self.sim_manager.pre_we()

    def test_post_we(self):
        self.sim_manager.post_we()
