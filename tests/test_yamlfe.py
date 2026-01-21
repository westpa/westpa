import pytest
import numpy as np
from numpy.testing import assert_array_equal

import westpa
import westpa.core.yamlcfg as ycf
from westpa.core.systems import WESTSystem
from westpa.core.binning import RectilinearBinMapper


class TESTSystem(WESTSystem):
    def initialize(self):
        self.pcoord_ndim = 1
        self.pcoord_dtype = np.float32
        self.pcoord_len = 5
        self.bin_mapper = RectilinearBinMapper([list(np.arange(0.0, 10.1, 0.1))])
        self.bin_target_counts = np.empty((self.bin_mapper.nbins,), np.int_)
        self.bin_target_counts[...] = 10
        self.test_variable_2 = "And I'm the second one"


# YAML Front end tests

# Implemented basic tests
#  - returns the correct system
#    given a system driver
#  - returns the correct system
#    given a yaml system
#  - returns the correct system
#    given both


# A class to test both paths at the same time
# if it works we assure we can load the driver
# AND overwrite it properly
class TestYAMLFrontEnd:
    def testYAMLFEDriver(self):
        '''
        Test method to ensure the YAML system generator works as
        advertised
        '''

        # First the objects that will be used for testing
        rc = westpa.rc
        yamlConf = ycf.YAMLConfig()
        # A sample dictionary from a yaml file
        test_dict = {
            "west": {
                "system": {
                    "driver": "tests.test_yamlfe.TESTSystem",
                    "system_options": {
                        "pcoord_ndim": 2,
                        "test_variable": "I'm a test variable",
                        "pcoord_len": 10,
                        "pcoord_dtype": np.float32,
                        "bin_target_counts": 10,
                        "bins": {"type": "RectilinearBinMapper", "boundaries": [[0.0, 0.5, 1.5, 2.5, 3.5, 'inf']]},
                    },
                }
            }
        }
        yamlConf._data = test_dict
        rc.config = yamlConf

        self.system = rc.new_system_driver()

        system = self.system
        # Assert we have the right options
        # This needs some more documentation and alerts for the assertions
        assert system.pcoord_ndim == 2
        assert system.test_variable == "I'm a test variable"
        # This one in particular checks if the bins are passed correctly
        assert (system.bin_mapper.boundaries == np.array([[0.0, 0.5, 1.5, 2.5, 3.5, 'inf']], dtype=np.float32)).all()
        assert system.pcoord_len == 10
        assert system.pcoord_dtype == np.float32
        ## These should be the same as the original
        assert system.test_variable_2 == "And I'm the second one"

    def testYAMLFEConfig(self):
        # First the objects that will be used for testing
        rc = westpa.rc
        yamlConf = ycf.YAMLConfig()
        # A sample dictionary from a yaml file
        test_dict = {
            "west": {
                "system": {
                    "system_options": {
                        "pcoord_ndim": 2,
                        "test_variable": "I'm a test variable",
                        "pcoord_len": 10,
                        "pcoord_dtype": np.float32,
                        "bin_target_counts": 10,
                        "bins": {"type": "RectilinearBinMapper", "boundaries": ["np.arange(0.0, 5.0, 0.5)"]},
                    }
                }
            }
        }
        yamlConf._data = test_dict
        rc.config = yamlConf

        self.system = rc.new_system_driver()

        system = self.system
        # Assert we have the right options
        # This needs some more documentation and alerts for the assertions
        assert system.pcoord_ndim == 2
        assert system.test_variable == "I'm a test variable"
        # This one in particular checks if the bins are passed correctly
        assert (system.bin_mapper.boundaries == np.arange(0.0, 5.0, 0.5)).all()
        assert system.pcoord_len == 10
        assert system.pcoord_dtype == np.float32

        # Test __setitem__() method of YAMLConfig()
        rc.config['west', 'propagation', 'max_total_iteration'] = 1000

        assert rc.config['west', 'propagation', 'max_total_iteration'] == 1000

    def testSystemDefaults(self):
        # First the objects that will be used for testing
        testSystem = WESTSystem()

        with pytest.raises(NotImplementedError):
            testSystem.new_region_set()

        # Test that the new pcoord array is of the correct shape
        test_zero = np.zeros((2, 1), np.float32)
        pcoord_array = testSystem.new_pcoord_array()
        assert testSystem.pcoord_len == 2
        assert_array_equal(test_zero, pcoord_array)

        testSystem.initialize()
        testSystem.prepare_run()
        testSystem.finalize_run()


class TestYAMLConfig:
    def test_dubious_config_entry(self):
        with pytest.warns(ycf.ConfigValueWarning):
            ycf.warn_dubious_config_entry('1', 1, expected_type=str)
            ycf.warn_dubious_config_entry('1', 1)

    def test_check_bool(self):
        with pytest.warns(ycf.ConfigValueWarning):
            ycf.check_bool(100, action='warn')
        with pytest.raises(ValueError):
            ycf.check_bool(100, action='raise')
