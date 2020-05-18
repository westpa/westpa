
from work_managers.threads import ThreadsWorkManager
from .tsupport import *

import nose.tools
from nose.tools import raises

class TestThreadsWorkManager(CommonWorkManagerTests,CommonParallelTests):
    def setUp(self):
        self.work_manager = ThreadsWorkManager()
        self.work_manager.startup()
        
    def tearDown(self):
        self.work_manager.shutdown()

class TestThreadsWorkManagerAux:
    def test_shutdown(self):
        work_manager = ThreadsWorkManager()
        work_manager.startup()
        work_manager.shutdown()
        for worker in work_manager.workers:
            assert not worker.is_alive() 

