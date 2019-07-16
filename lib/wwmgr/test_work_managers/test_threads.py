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

