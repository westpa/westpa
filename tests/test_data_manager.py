import argparse
import os

import numpy as np
import westpa

""" IDEA: instantiate data manager and check that everything is set up correctly
    question: what state is the data manager expected to be in?
    also: look at data manger file to see other basic functions to call to check
    that set up is working correctly
    * and need to be careful since apparently DataManager state is persisting after
    sim manager test is run? <-- need to check on this and make sure we don't have similar issues here
"""


class TestDataManager:

    def setup(self):
        """ Literally just writing what Josh and John did to set up rc stuff
         in their sim manager test"""
        parser = argparse.ArgumentParser()
        westpa.rc.add_args(parser)
        here = os.path.dirname(__file__)
        # Set SIM_ROOT to fixtures folder with west.cfg for odld simulation
        os.environ['WEST_SIM_ROOT'] = os.path.join(here, 'fixtures', 'odld')

        config_file_name = os.path.join(here, 'fixtures', 'odld', 'west.cfg')
        args = parser.parse_args(['-r={}'.format(config_file_name)])
        westpa.rc.process_args(args)

        self.data_manager = westpa.rc.get_data_manager()

    def teardown(self):
        # Also modeled after sim_manager test teardown
        westpa.rc._data_manager = None
        westpa.rc._system = None
        del os.environ['WEST_SIM_ROOT']

    def test_data_manager(self):
        assert(len(self.data_manager.get_basis_states()) == 0)

        self.get_iter_summary()

    def test_target_states(self):
        self.data_manager.get_target_states()

    def test_inital_states(self):
        assert(len(self.data_manager.get_initial_states()) == 0)
