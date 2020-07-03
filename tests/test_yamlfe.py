import numpy as np

import westpa
import westpa.core.yamlcfg as ycf
from westpa.core.binning import RectilinearBinMapper


class TESTSystem(ycf.YAMLSystem):
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
