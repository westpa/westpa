# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.


# Need to sort out exactly what I need to import here
import westpa, numpy, copy
import westpa.yamlcfg as ycf

from westpa.binning import RectilinearBinMapper

import logging
log = logging.getLogger(__name__)


class TESTSystem(ycf.YAMLSystem):
    def initialize(self):
        self.pcoord_ndim = 1
        self.pcoord_dtype = numpy.float32
        self.pcoord_len = 5
        self.bin_mapper = RectilinearBinMapper([ list(numpy.arange(0.0, 10.1, 0.1)) ] )
        self.bin_target_counts = numpy.empty((self.bin_mapper.nbins,), numpy.int_)
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
        test_dict = {"west":{
                       "system": {
                         "driver" :"testyamlfe.TESTSystem",
                         "system_options": { 
                           "pcoord_ndim":2, "test_variable": "I'm a test variable",
                           "pcoord_len":10, "pcoord_dtype": numpy.float32,
                           "bin_target_counts": 10,
                             "bins": {
                               "type":"RectilinearBinMapper", 
                                 "boundaries": [[0.0, 0.5, 1.5, 2.5, 3.5, 'inf']]
                                 }}}}}
        yamlConf._data = test_dict
        rc.config = yamlConf

        self.system = rc.new_system_driver()
       
        system = self.system
        # Assert we have the right options
        # This needs some more documentation and alerts for the assertions
        assert system.pcoord_ndim == 2
        assert system.test_variable == "I'm a test variable"
        # This one in particular checks if the bins are passed correctly
        assert (system.bin_mapper.boundaries == \
               numpy.array([[0.0, 0.5, 1.5, 2.5, 3.5, 'inf']], dtype=numpy.float32)).all()
        assert system.pcoord_len == 10
        assert system.pcoord_dtype == numpy.float32
        ## These should be the same as the original
        assert system.test_variable_2 == "And I'm the second one"

    def testYAMLFEConfig(self):
        # First the objects that will be used for testing
        rc = westpa.rc
        yamlConf = ycf.YAMLConfig()    
        # A sample dictionary from a yaml file
        test_dict = {"west":{
                       "system": {
                         "system_options": { 
                           "pcoord_ndim":2, "test_variable": "I'm a test variable",
                           "pcoord_len":10, "pcoord_dtype": numpy.float32,
                           "bin_target_counts": 10,
                             "bins": {
                               "type":"RectilinearBinMapper", 
                                 "boundaries": ["numpy.arange(0.0, 5.0, 0.5)"]
                                 }}}}}
        yamlConf._data = test_dict
        rc.config = yamlConf

        self.system = rc.new_system_driver()

        system = self.system
        # Assert we have the right options
        # This needs some more documentation and alerts for the assertions
        assert system.pcoord_ndim == 2
        assert system.test_variable == "I'm a test variable"
        # This one in particular checks if the bins are passed correctly
        assert (system.bin_mapper.boundaries == \
               numpy.arange(0.0, 5.0, 0.5)).all()
        assert system.pcoord_len == 10
        assert system.pcoord_dtype == numpy.float32
