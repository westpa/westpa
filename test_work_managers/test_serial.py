'''
Created on Apr 26, 2012

@author: mzwier

Tests work manager base class via serial interface
'''

from work_managers.serial import SerialWorkManager
from tsupport import *

import nose.tools
from nose.tools import raises

class TestSerialWorkManager(CommonWorkManagerTests):
    def setUp(self):
        self.work_manager = SerialWorkManager()