import unittest

from westpa.work_managers.mpi import MPIWorkManager

# from .tsupport import CommonWorkManagerTests, CommonParallelTests


# class TestMPIWorkManager(CommonWorkManagerTests, CommonParallelTests):
class TestMPIWorkManager(unittest.TestCase):
    def setUp(self):
        self.work_manager = MPIWorkManager()
        self.work_manager.startup()

    def tearDown(self):
        self.work_manager.shutdown()

    def test_null(self):
        assert 1 == 1
