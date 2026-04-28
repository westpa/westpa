import unittest

from westpa.work_managers.threads import ThreadsWorkManager
from .tsupport import CommonWorkManagerTests, CommonParallelTests


class TestThreadsWorkManager(unittest.TestCase, CommonWorkManagerTests, CommonParallelTests):
    @classmethod
    def setUpClass(cls):
        cls.work_manager = ThreadsWorkManager(n_workers=5)
        cls.work_manager.startup()

    @classmethod
    def tearDownClass(cls):
        cls.work_manager.shutdown()


class TestThreadsWorkManagerAux:
    def test_shutdown(self):
        work_manager = ThreadsWorkManager(n_workers=5)
        work_manager.startup()
        work_manager.shutdown()
        for worker in work_manager.workers:
            assert not worker.is_alive()
