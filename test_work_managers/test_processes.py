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
        
    @nose.tools.timed(2)                    
    def test_worker_ids(self):
        work_manager = ProcessWorkManager()
        with work_manager:
            futures = work_manager.submit_many([(get_process_index, (), {})] * work_manager.n_workers)
            work_manager.wait_all(futures)
            results = set(future.get_result() for future in futures)
            assert results == set(str(n) for n in xrange(work_manager.n_workers)), results
        
        