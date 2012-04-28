'''
Created on Apr 26, 2012

@author: mzwier

Tests work manager base class via serial interface
'''

from work_managers.serial import SerialWorkManager
from tsupport import *

import nose.tools
from nose.tools import raises

    
class TestSerialWorkManager:
    def test_submit(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            
    def test_submit_many(self):
        with SerialWorkManager() as work_manager:
            futures = work_manager.submit([(will_succeed,(),{}) for i in xrange(5)])
            
    def test_as_completed(self):
        with SerialWorkManager() as work_manager:
            input = set(xrange(5))
            futures = [work_manager.submit(identity, i) for i in xrange(5)]
            output = set(future.get_result() for future in work_manager.as_completed(futures))
            assert input == output
        
    def test_wait_any(self):
        with SerialWorkManager() as work_manager:
            input = set(xrange(5))
            futures = [work_manager.submit(identity, i) for i in xrange(5)]
            output = work_manager.wait_any(futures).get_result()
            assert output in input
        
    def test_wait_all(self):
        with SerialWorkManager() as work_manager:
            input = set(xrange(5))
            futures = [work_manager.submit(identity, i) for i in xrange(5)]
            output = set(future.get_result() for future in work_manager.wait_all(futures))
            assert input == output
        

class TestWMFuture:        
    def test_result(self):
        with SerialWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            assert future.get_result() is True    
    
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
        
   
