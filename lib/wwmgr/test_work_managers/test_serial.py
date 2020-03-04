
from work_managers.serial import SerialWorkManager
from .tsupport import *

import nose.tools
from nose.tools import raises

class TestSerialWorkManager(CommonWorkManagerTests):
    def setUp(self):
        self.work_manager = SerialWorkManager()