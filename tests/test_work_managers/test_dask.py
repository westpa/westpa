import pytest

from westpa.work_managers import DaskWorkManager


@pytest.fixture(scope='module')
def work_manager():
    return DaskWorkManager()


def test_submit(work_manager):
    future = work_manager.submit(sum, args=[(1, 2)])
    assert not future.done
    assert future.result == 3
    assert future.done


def test_as_completed(work_manager):
    futures = [work_manager.submit(str, args=[index]) for index in range(10)]
    for future in work_manager.as_completed(futures):
        print(future)
        assert future.done
        assert future.result == str(futures.index(future))


def test_submit_as_completed(work_manager):
    tasks = [(str, [index], None) for index in range(10)]
    for future in work_manager.submit_as_completed(iter(tasks), queue_size=4):
        assert future.done
        assert isinstance(future.result, str)


def test_wait_any(work_manager):
    futures = [work_manager.submit(str, args=[index]) for index in range(10)]
    future = work_manager.wait_any(futures)
    assert future.done
    assert future.result == str(futures.index(future))


def test_exception(work_manager):
    future = work_manager.submit(int, args=['a'])
    with pytest.raises(ValueError):
        _ = future.result
    assert isinstance(future.exception, ValueError)
    assert future.done
