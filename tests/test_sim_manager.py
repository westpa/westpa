import argparse
import os

import numpy as np

import westpa
from westpa.core.binning.assign import RectilinearBinMapper


class TestSimManager:
    def setup(self):

        parser = argparse.ArgumentParser()
        westpa.rc.add_args(parser)

        here = os.path.dirname(__file__)
        os.environ['WEST_SIM_ROOT'] = os.path.join(here, 'fixtures', 'odld')

        config_file_name = os.path.join(here, 'fixtures', 'odld', 'west.cfg')
        args = parser.parse_args(['-r={}'.format(config_file_name)])
        westpa.rc.process_args(args)
        self.sim_manager = westpa.rc.get_sim_manager()

    def teardown(self):
        westpa.rc._sim_manager = None
        westpa.rc._system = None
        del os.environ['WEST_SIM_ROOT']

    def test_sim_manager(self):
        assert self.sim_manager.n_propagated == 0
        assert len(self.sim_manager._callback_table) == 0

    def dummy_callback_one(self):
        system = self.sim_manager.system

        bounds = [0.0, 1.0, 2.0, 3.0]
        system.bin_mapper = RectilinearBinMapper([bounds])

    def dummy_callback_two(self):
        system = self.sim_manager.system

        bounds = [0.0, 1.0, 2.0, 5.0]
        system.bin_mapper = RectilinearBinMapper([bounds])

    def test_register_callback(self):
        hook = self.sim_manager.prepare_new_iteration

        self.sim_manager.register_callback(hook, self.dummy_callback_one, 3)
        self.sim_manager.register_callback(hook, self.dummy_callback_two, 0)
        assert hook in self.sim_manager._callback_table

        callbacks = self.sim_manager._callback_table.get(hook, [])
        assert (3, self.dummy_callback_one.__name__, self.dummy_callback_one) in callbacks
        assert (0, self.dummy_callback_two.__name__, self.dummy_callback_two) in callbacks

    def test_invoke_callback(self):
        hook = self.sim_manager.prepare_new_iteration

        self.sim_manager.register_callback(hook, self.dummy_callback_one, 3)
        self.sim_manager.register_callback(hook, self.dummy_callback_two, 0)

        self.sim_manager.invoke_callbacks(hook)

        system = self.sim_manager.system
        assert np.all(system.bin_mapper.boundaries == np.array([0.0, 1.0, 2.0, 3.0]))
