import unittest

from westpa.work_managers.serial import SerialWorkManager
from .tsupport import CommonWorkManagerTests


class TestSerialWorkManager(unittest.TestCase, CommonWorkManagerTests):
    def setUp(self):
        self.work_manager = SerialWorkManager()
