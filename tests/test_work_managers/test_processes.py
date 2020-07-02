import os
import signal
import unittest

from westpa.work_managers.processes import ProcessWorkManager
from .tsupport import CommonParallelTests, CommonWorkManagerTests
from .tsupport import will_busyhang, will_busyhang_uninterruptible, get_process_index

import nose.tools
from nose.tools import raises


class TestProcessWorkManager(unittest.TestCase, CommonParallelTests, CommonWorkManagerTests):
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
        for i in range(5):
            work_manager.submit(will_busyhang)
        work_manager.shutdown()
        for worker in work_manager.workers:
            assert not worker.is_alive()

    @nose.tools.timed(2)
    def test_hang_shutdown_ignoring_sigint(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()
        for i in range(5):
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
        for i in range(5):
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
            assert results == set(str(n) for n in range(work_manager.n_workers)), results
