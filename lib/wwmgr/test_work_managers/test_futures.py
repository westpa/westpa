
from work_managers.serial import SerialWorkManager
from nose.tools import assert_raises #@UnresolvedImport
from .tsupport import *

class TestWMFuture:        
    def test_result(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            assert future.get_result() is True
            
    def test_discarded_result(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            assert future.get_result(discard=True) is True
            assert_raises(AttributeError, getattr, future, '_result')
    
    @raises(ExceptionForTest)
    def test_exception_raise(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_fail)
            future.get_result()
        
    def test_exception_retrieve(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_fail)
            exc = future.get_exception()
            assert exc.args[0] == 'failed as expected'
        
    def test_callback(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            def cbfn(future):
                assert future.get_result() is True
            future._add_callback(cbfn)
        
    def test_success_wait(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            future.wait()
    
    def test_exception_wait(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_fail)
            future.wait()
        
    def test_is_done(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            future.wait()
            assert future.done
        
