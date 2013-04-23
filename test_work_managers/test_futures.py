# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

from work_managers.serial import SerialWorkManager
from nose.tools import assert_raises #@UnresolvedImport
from tsupport import *

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
        
