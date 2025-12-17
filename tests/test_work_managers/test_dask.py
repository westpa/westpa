import pytest

from westpa.work_managers import DaskWorkManager


def test_submit():
    work_manager = DaskWorkManager()
    future = work_manager.submit(sum, args=[(1, 2)])
    assert not future.done
    assert future.result == 3
    assert future.done


def test_as_completed():
    work_manager = DaskWorkManager()
    futures = [work_manager.submit(str, args=[index]) for index in range(10)]
    for future in work_manager.as_completed(futures):
        assert future.result == str(futures.index(future))


def test_exception():
    work_manager = DaskWorkManager()
    future = work_manager.submit(int, args=['a'])
    with pytest.raises(ValueError):
        _ = future.result
    assert isinstance(future.exception, ValueError)
    assert future.done
