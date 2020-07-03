import argparse
import os
import unittest

import westpa.work_managers.environment
from westpa.work_managers.environment import make_work_manager, add_wm_args, process_wm_args
from westpa.work_managers import SerialWorkManager, ThreadsWorkManager, ProcessWorkManager, ZMQWorkManager

from .tsupport import will_succeed


class TestInstantiations(unittest.TestCase):
    '''Test to see that the environment system at least selects the proper work manager and sets the
    number of workers appropriately'''

    sanitize_vars = ('WM_WORK_MANAGER', 'WM_N_WORKERS')

    def setUp(self):
        for varname in self.sanitize_vars:
            assert varname not in os.environ
        westpa.work_managers.environment.default_env.args = None

    def tearDown(self):
        for varname in self.sanitize_vars:
            os.environ.pop(varname, None)
        westpa.work_managers.environment.default_env.args = None

    def testArgs(self):
        parser = argparse.ArgumentParser()
        add_wm_args(parser)
        args = '--work-manager=threads --n-workers=3'.split()
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
        args = '--work-manager=threads'.split()
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
        with work_manager:
            # Need to send enough work to start sufficient workers
            for _ in range(2):
                future = work_manager.submit(will_succeed)
                future.get_result()

            assert work_manager.n_workers == 3
