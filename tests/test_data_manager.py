import unittest
import argparse
import os
import tempfile
import westpa


class TestDataManager(unittest.TestCase):
    test_name = 'DATA_MANAGER'

    def __init__(self, methodName):
        super().__init__(methodName='test_data_manager')

    def setUp(self):
        parser = argparse.ArgumentParser()
        westpa.rc.add_args(parser)

        here = os.path.dirname(__file__)
        # Set SIM_ROOT to fixtures folder with west.cfg for odld simulation
        os.environ['WEST_SIM_ROOT'] = os.path.join(here, 'fixtures', 'odld')

        config_file_name = os.path.join(here, 'fixtures', 'odld', 'west.cfg')
        args = parser.parse_args(['-r={}'.format(config_file_name), "--verbose"])
        westpa.rc.process_args(args)

        self.sim_manager = westpa.rc.get_sim_manager()
        self.test_dir = tempfile.mkdtemp()
        self.hdf5 = os.path.join(self.test_dir, "west.h5")
        # data = westpa.rc.get_sim_manager().we_driver.rc.get_data_manager()
        # data.we_h5filename = self.hdf5

        self.data_manager = westpa.rc.new_data_manager()
        self.data_manager.we_h5filename = self.hdf5
        self.assertEqual(self.data_manager.h5_access_mode, 'r+')
        """ 1. westpa.rc.get_data_manager() is executed and calls the new_data_manager function.
            2. The new_data_manager instantiates a WESTDataManager object, thus calling the
               WESTDataManager.__init__() constructor. In __init__, process_config() is executed.
            3. process_config() will read in the fixtures/odld/west.cfg.
               That's all that the setup does in it's current state.
            4. The test checks that the defaults that are set in __init__ are
               indeed the defaults and that the data-manager relevant parts specified
               in west.cfg are indeed updated accordingly."""

    def tearDown(self):
        del os.environ['WEST_SIM_ROOT']
        westpa.rc = westpa.core._rc.WESTRC()

    def test_data_manager(self):
        assert self.data_manager.h5_access_mode == 'r+'
        assert self.data_manager.closed
        assert os.path.basename(self.data_manager.we_h5filename) == 'west.h5'
        assert self.data_manager.aux_compression_threshold == 16384
        assert len(self.data_manager.dataset_options) == 2
