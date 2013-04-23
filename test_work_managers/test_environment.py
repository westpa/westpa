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

import os, argparse
import work_managers.environment
from work_managers.environment import make_work_manager, add_wm_args, process_wm_args
from work_managers import SerialWorkManager, ThreadsWorkManager, ProcessWorkManager, ZMQWorkManager
from tsupport import will_succeed

class TestInstantiations:
    '''Test to see that the environment system at least selects the proper work manager and sets the
    number of workers appropriately'''
    
    sanitize_vars = ('WM_WORK_MANAGER', 'WM_N_WORKERS')
    
    def setUp(self):
        for varname in self.sanitize_vars:
            assert varname not in os.environ
        work_managers.environment.default_env.args = None            
    
    def tearDown(self):
        for varname in self.sanitize_vars:
            os.environ.pop(varname, None)
        work_managers.environment.default_env.args = None
            
    def testArgs(self):
        parser = argparse.ArgumentParser()
        add_wm_args(parser)
        args='--work-manager=threads --n-workers=3'.split()
        args = parser.parse_args(args)
        process_wm_args(args)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ThreadsWorkManager)
        assert work_manager.n_workers == 3
        
    def testArgFallthrough(self):
        # this test specifies the work manager on the command line, and the worker count in the environment
        # this simply tests whether we look to the environment in the case of a missing command line argument
        os.environ['WM_N_WORKERS'] = str(3)
        parser = argparse.ArgumentParser()
        add_wm_args(parser)
        args='--work-manager=threads'.split()
        args = parser.parse_args(args)
        process_wm_args(args)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ThreadsWorkManager)
        assert work_manager.n_workers == 3
        
    def testSerial(self):
        os.environ['WM_WORK_MANAGER'] = 'serial'
        work_manager = make_work_manager()
        assert isinstance(work_manager, SerialWorkManager)
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        
    def testThreads(self):
        os.environ['WM_WORK_MANAGER'] = 'threads'
        os.environ['WM_N_WORKERS'] = str(3)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ThreadsWorkManager)
        assert work_manager.n_workers == 3
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        
    def testProcesses(self):
        os.environ['WM_WORK_MANAGER'] = 'processes'
        os.environ['WM_N_WORKERS'] = str(3)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ProcessWorkManager)
        assert work_manager.n_workers == 3
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        
    def testZeroMQ(self):
        os.environ['WM_WORK_MANAGER'] = 'zmq'
        os.environ['WM_N_WORKERS'] = str(3)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ZMQWorkManager)
        assert work_manager.internal_client.n_workers == 3
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        