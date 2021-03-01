import argparse
import os
from unittest import TestCase
from unittest.mock import MagicMock
import tempfile
import pytest

import numpy as np

import westpa
from westpa.core.binning.assign import RectilinearBinMapper
from westpa.core.segment import Segment
from westpa.core.states import BasisState
from westpa.core.sim_manager import PropagationError


class TestSimManager(TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        westpa.rc.add_args(parser)

        here = os.path.dirname(__file__)
        os.environ['WEST_SIM_ROOT'] = os.path.join(here, 'fixtures', 'odld')

        config_file_name = os.path.join(here, 'fixtures', 'odld', 'west.cfg')
        args = parser.parse_args(['-r={}'.format(config_file_name)])
        westpa.rc.process_args(args)
        self.sim_manager = westpa.rc.get_sim_manager()
        self.test_dir = tempfile.mkdtemp()
        self.hdf5 = os.path.join(self.test_dir, "west.h5")
        self.basis_states = [BasisState(label="label", probability=1.0)]
        self.segments = [self.segment(0.0, 1.5, weight=0.125) for _i in range(4)] + [
            self.segment(1.5, 0.5, weight=0.125) for _i in range(4)
        ]
        self.sim_manager.we_driver.new_iteration()
        self.sim_manager.we_driver.assign(self.segments)
        self.sim_manager.we_driver.construct_next()
        self.sim_manager.segments = {segment.seg_id: segment for segment in self.segments}
        self.sim_manager.incomplete_segments = self.sim_manager.segments
        self.sim_manager.current_iter_istates = self.sim_manager.segments
        self.sim_manager.completed_segments = self.sim_manager.segments
        self.sim_manager.report_bin_statistics = MagicMock(return_value=True)

        data = self.sim_manager.we_driver.rc.get_data_manager()
        data.we_h5filename = self.hdf5
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
        data.get_segments = MagicMock(return_value=self.segments)
        self.sim_manager.we_driver.rc.get_data_manager = MagicMock(return_value=data)
        self.sim_manager.n_iter = n_iter

    def tearDown(self):
        westpa.rc._sim_manager = None
        westpa.rc._system = None
        westpa.rc._data_manager = None
        del os.environ['WEST_SIM_ROOT']

    def dummy_callback_one(self):
        system = self.sim_manager.system
        bounds = [0.0, 1.0, 2.0, 3.0]
        system.bin_mapper = RectilinearBinMapper([bounds])

    def dummy_callback_two(self):
        system = self.sim_manager.system
        bounds = [0.0, 1.0, 2.0, 5.0]
        system.bin_mapper = RectilinearBinMapper([bounds])

    def segment(self, init_pcoord, final_pcoord, weight=1.0):
        segment = Segment(n_iter=1, seg_id=1123, pcoord=self.sim_manager.system.new_pcoord_array(), weight=weight)
        segment.pcoord[0] = init_pcoord
        segment.pcoord[1] = final_pcoord
        segment.parent_id = 1
        return segment

    def test_sim_manager(self):
        self.assertEquals(self.sim_manager.n_propagated, 0)
        self.assertEquals(len(self.sim_manager._callback_table), 0)

    def test_register_callback(self):
        hook = self.sim_manager.prepare_new_iteration

        self.sim_manager.register_callback(hook, self.dummy_callback_one, 3)
        self.sim_manager.register_callback(hook, self.dummy_callback_two, 0)
        self.assertTrue(hook in self.sim_manager._callback_table)

        callbacks = self.sim_manager._callback_table.get(hook, [])
        self.assertTrue((3, self.dummy_callback_one.__name__, self.dummy_callback_one) in callbacks)  # noqa
        self.assertTrue((0, self.dummy_callback_two.__name__, self.dummy_callback_two) in callbacks)  # noqa

    def test_invoke_callback(self):
        hook = self.sim_manager.prepare_new_iteration

        self.sim_manager.register_callback(hook, self.dummy_callback_one, 3)
        self.sim_manager.register_callback(hook, self.dummy_callback_two, 0)

        self.sim_manager.invoke_callbacks(hook)

        system = self.sim_manager.system
        self.assertTrue(np.all(system.bin_mapper.boundaries == np.array([0.0, 1.0, 2.0, 3.0])))  # noqa

    def test_process_config(self):
        self.sim_manager.process_config()
        self.assertTrue(self.sim_manager.do_gen_istates)
        self.assertEquals(self.sim_manager.propagator_block_size, 10000)
        self.assertFalse(self.sim_manager.save_transition_matrices)
        self.assertEquals(self.sim_manager.max_run_walltime, 10800)
        self.assertEquals(self.sim_manager.max_total_iterations, 100)

    def test_load_plugins(self):
        self.sim_manager.load_plugins()

    def test_report_bin_statistics(self):
        self.sim_manager.report_bin_statistics([0.0, 1.0, 2.0, 5.0])

    def test_get_bstate_pcoords(self):
        self.sim_manager.get_bstate_pcoords(self.basis_states)

    def test_report_basis_states(self):
        self.sim_manager.report_basis_states(self.basis_states)

    def test_report_target_states(self):
        self.sim_manager.report_target_states(self.basis_states)

    @pytest.mark.skip(reason="Cannot currently test WESimManager.initialize_simulation()")
    def test_initialize_simulation(self):
        # TODO: determine how to test self.simulation_manager.initialize_simulation()
        pass

    def test_prepare_iteration(self):
        self.sim_manager.prepare_new_iteration()
        self.sim_manager.prepare_iteration()

    def test_finalize_iteration(self):
        self.sim_manager.finalize_iteration()

    def test_get_istate_futures(self):
        self.sim_manager.get_istate_futures()

    def test_propagate(self):
        westpa.core.states.pare_basis_initial_states = MagicMock(return_value=([], []))
        self.sim_manager.propagate

    def test_save_bin_data(self):
        self.sim_manager.save_bin_data()

    def test_check_propagation(self):
        self.assertRaises(PropagationError, self.sim_manager.check_propagation)

    def test_run_we(self):
        self.sim_manager.run_we()

    def test_run(self):
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
