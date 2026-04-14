import dask.distributed as distributed
import pytest

from westpa.work_managers import DaskWorkManager

NUM_TASKS = 5
QUEUE_SIZE = 2


class TestDaskWorkManager:

    @pytest.fixture(scope='class')
    def work_manager(self):
        # Use a faster, lighter cluster for testing
        cluster = distributed.LocalCluster(n_workers=3, threads_per_worker=1, memory_limit='1GB')
        client = distributed.Client(cluster)
        with DaskWorkManager(client) as work_manager:
            yield work_manager

        # Cleanup done during fixture teardown
        client.close()
        cluster.close()

    def test_submit(self, work_manager):
        future = work_manager.submit(sum, args=[(1, 2)])
        assert not future.done
        assert future.result == 3
        assert future.done

    def test_as_completed(self, work_manager):
        futures = [work_manager.submit(str, args=[index]) for index in range(NUM_TASKS)]
        for future in work_manager.as_completed(futures):
            print(future)
            assert future.done
            assert future.result == str(futures.index(future))

    def test_submit_as_completed(self, work_manager):
        tasks = [(str, [index], None) for index in range(NUM_TASKS)]
        for future in work_manager.submit_as_completed(iter(tasks), queue_size=QUEUE_SIZE):
            assert future.done
            assert isinstance(future.result, str)

    def test_wait_any(self, work_manager):
        futures = [work_manager.submit(str, args=[index]) for index in range(NUM_TASKS)]
        future = work_manager.wait_any(futures)
        assert future.done
        assert future.result == str(futures.index(future))

    def test_exception(self, work_manager):
        future = work_manager.submit(int, args=['a'])
        with pytest.raises(ValueError):
            _ = future.result
        assert isinstance(future.exception, ValueError)
        assert future.done
