import argparse
import os

import numpy as np
import westpa
import logging


class TestDataManager:

    def setup(self):
        parser = argparse.ArgumentParser()
        westpa.rc.add_args(parser)

        here = os.path.dirname(__file__)
        # Set SIM_ROOT to fixtures folder with west.cfg for odld simulation
        os.environ['WEST_SIM_ROOT'] = os.path.join(here, 'fixtures', 'odld')

        config_file_name = os.path.join(here, 'fixtures', 'odld', 'west.cfg')
        args = parser.parse_args(['-r={}'.format(config_file_name), "--verbose"])
        westpa.rc.process_args(args)

        self.data_manager = westpa.rc.get_data_manager()

        """ 1. westpa.rc.get_data_manager() is executed and calls the new_data_manager function.
            2. The new_data_manager instantiates a WESTDataManager object, thus calling the
               WESTDataManager.__init__() constructor. In __init__, process_config() is executed.
            3. process_config() will read in the fixtures/odld/west.cfg.
               That's all that the setup does in it's current state.
            4. The test checks that the defaults that are set in __init__ are
               indeed the defaults and that the data-manager relevant parts specified
               in west.cfg are indeed updated accordingly."""

    def teardown(self):
        westpa.rc._data_manager = None
        westpa.rc._system = None
        del os.environ['WEST_SIM_ROOT']

    def test_data_manager(self):
        assert self.data_manager.h5_access_mode == 'r+'
        assert self.data_manager.we_h5file == None
        assert os.path.basename(self.data_manager.we_h5filename) == 'west.h5'
        assert self.data_manager.aux_compression_threshold == 16384
        assert len(self.data_manager.dataset_options) == 2
