"""import os, signal

from work_managers.processes import ProcessWorkManager
from tsupport import *

import nose.tools
from nose.tools import raises

class TestProcessWorkManager:
    def test_submit(self):
        with ProcessWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)

    def test_submit_many(self):
        with ProcessWorkManager() as work_manager:
            futures = work_manager.submit([(will_succeed,(),{}) for i in xrange(5)])
            
    @nose.tools.timed(10)
    def test_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.startup()
        work_manager.shutdown()
        for worker in work_manager.workers:
            assert not worker.is_alive()

    @nose.tools.timed(10)            
    def test_hang_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 1
        work_manager.startup()
        for i in xrange(5):
            work_manager.submit(will_busyhang)
        work_manager.shutdown() 
        for worker in work_manager.workers:
            assert not worker.is_alive()

    @nose.tools.timed(10)            
    def test_hang_shutdown_ignoring_sigint(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 1
        work_manager.startup()
        for i in xrange(5):
            work_manager.submit(will_busyhang_uninterruptible)
        work_manager.shutdown() 
        for worker in work_manager.workers:
            assert not worker.is_alive()
            
    @nose.tools.timed(10)
    @raises(KeyboardInterrupt)
    def test_sigint_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 1
        work_manager.startup()
        for i in xrange(5):
            work_manager.submit(will_busyhang)
    
        try:
            os.kill(os.getpid(), signal.SIGINT)
        except KeyboardInterrupt:
            for worker in work_manager.workers:
                assert not worker.is_alive()
            raise
                    
    def test_as_completed(self):
        with ProcessWorkManager() as work_manager:
            input = set(xrange(5))
            futures = [work_manager.submit(identity, i) for i in xrange(5)]
            output = set(future.get_result() for future in work_manager.as_completed(futures))
            assert input == output
        
    def test_wait_any(self):
        with ProcessWorkManager() as work_manager:
            input = set(xrange(5))
            futures = [work_manager.submit(identity, i) for i in xrange(5)]
            output = work_manager.wait_any(futures).get_result()
            assert output in input
        
    def test_wait_all(self):
        with ProcessWorkManager() as work_manager:
            input = set(xrange(5))
            futures = [work_manager.submit(identity, i) for i in xrange(5)]
            output = set(future.get_result() for future in work_manager.wait_all(futures))
            assert input == output

    def test_result(self):
        with ProcessWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            assert future.get_result() is True    
    
    @raises(ExceptionForTest)
    def test_exception_raise(self):
        with ProcessWorkManager() as work_manager:
            future = work_manager.submit(will_fail)
            future.get_result()
        
    def test_exception_retrieve(self):
        with ProcessWorkManager() as work_manager:
            future = work_manager.submit(will_fail)
            exc = future.get_exception()
            assert exc.args[0] == 'failed as expected'
                
    def test_success_wait(self):
        with ProcessWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            future.wait()
    
    def test_exception_wait(self):
        with ProcessWorkManager() as work_manager:
            future = work_manager.submit(will_fail)
            future.wait()
        
    def test_is_done(self):
        with ProcessWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            future.wait()
            assert future.done
"""