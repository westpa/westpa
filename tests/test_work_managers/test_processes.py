import os
import signal
import unittest
import pytest

from westpa.work_managers.processes import ProcessWorkManager
from .tsupport import CommonParallelTests, CommonWorkManagerTests
from .tsupport import will_busyhang, will_busyhang_uninterruptible, get_process_index


class TestProcessWorkManager(unittest.TestCase, CommonParallelTests, CommonWorkManagerTests):
    def setUp(self):
        self.work_manager = ProcessWorkManager()
        self.work_manager.startup()

    def tearDown(self):
        self.work_manager.shutdown()


class TestProcessWorkManagerAux:
    @pytest.mark.timeout(5)
    def test_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.startup()
        work_manager.shutdown()
        for worker in work_manager.workers:
            try:
                assert not worker.is_alive()
            except ValueError:
                pass  # probably closed already

    @pytest.mark.timeout(5)
    def test_hang_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()
        for _ in range(5):
            work_manager.submit(will_busyhang)
        work_manager.shutdown()
        for worker in work_manager.workers:
            try:
                assert not worker.is_alive()
            except ValueError:
                pass  # probably closed already

    @pytest.mark.timeout(5)
    def test_hang_shutdown_ignoring_sigint(self):
        work_manager = ProcessWorkManager()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()
        for _ in range(5):
            work_manager.submit(will_busyhang_uninterruptible)
        work_manager.shutdown()
        for worker in work_manager.workers:
            try:
                assert not worker.is_alive()
            except ValueError:
                pass  # probably closed already

    @pytest.mark.timeout(5)
    def test_sigint_shutdown(self):
        work_manager = ProcessWorkManager()
        work_manager.install_sigint_handler()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()
        for _ in range(5):
            work_manager.submit(will_busyhang)

        with pytest.raises(KeyboardInterrupt):
            try:
                os.kill(os.getpid(), signal.SIGINT)
            except KeyboardInterrupt:
                for worker in work_manager.workers:
                    try:
                        assert not worker.is_alive()
                    except ValueError:
                        pass  # probably closed already
                raise

    @pytest.mark.timeout(5)
    def test_worker_close_fail(self, monkeypatch):
        work_manager = ProcessWorkManager()
        work_manager.install_sigint_handler()
        work_manager.shutdown_timeout = 0.1
        work_manager.startup()

        work_manager.submit(will_busyhang_uninterruptible)
        worker = work_manager.workers[0]

        with monkeypatch.context() as m:
            m.setattr(worker, 'close', lambda: exec('raise(ValueError)'))
            m.setattr(worker, 'is_alive', lambda: True)
            m.setattr(work_manager, '_empty_queues', lambda: 0)
            work_manager.shutdown()

        # Clean up
        work_manager.shutdown()

    @pytest.mark.timeout(5)
    def test_worker_ids(self):
        work_manager = ProcessWorkManager()
        with work_manager:
            futures = work_manager.submit_many([(get_process_index, (), {})] * work_manager.n_workers)
            work_manager.wait_all(futures)
            results = set(future.get_result() for future in futures)
            assert results == set(str(n) for n in range(work_manager.n_workers)), results
