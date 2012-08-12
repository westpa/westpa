import os, signal

from work_managers.processes import ProcessWorkManager
from tsupport import *

import nose.tools
from nose.tools import raises

class TestProcessWorkManager(CommonParallelTests,CommonWorkManagerTests):
    def setUp(self):
        self.work_manager = ProcessWorkManager()
        self.work_manager.startup()
    def tearDown(self):
        self.work_manager.shutdown()

class TestProcessWorkManagerAux:            
    @nose.tools.timed(2)
    def test_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.startup()
        work_manager.shutdown()
        for worker in work_manager.workers:
            assert not worker.is_alive()

    @nose.tools.timed(2)            
    def test_hang_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()
        for i in xrange(5):
            work_manager.submit(will_busyhang)
        work_manager.shutdown() 
        for worker in work_manager.workers:
            assert not worker.is_alive()

    @nose.tools.timed(2)            
    def test_hang_shutdown_ignoring_sigint(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()
        for i in xrange(5):
            work_manager.submit(will_busyhang_uninterruptible)
        work_manager.shutdown() 
        for worker in work_manager.workers:
            assert not worker.is_alive()
            
    @nose.tools.timed(2)
    @raises(KeyboardInterrupt)
    def test_sigint_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.install_sigint_handler()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()
        for i in xrange(5):
            work_manager.submit(will_busyhang)
    
        try:
            os.kill(os.getpid(), signal.SIGINT)
        except KeyboardInterrupt:
            for worker in work_manager.workers:
                assert not worker.is_alive()
            raise
                    
