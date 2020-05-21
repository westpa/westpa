from westpa.work_managers.serial import SerialWorkManager
from .tsupport import CommonWorkManagerTests


class TestSerialWorkManager(CommonWorkManagerTests):
    def setUp(self):
        self.work_manager = SerialWorkManager()
