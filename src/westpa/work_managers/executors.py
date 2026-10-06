"""Serial implementation of the :mod:`concurrent.futures` executor interface."""

from concurrent.futures import Executor, Future

__all__ = ['SerialExecutor']


def _iter_outcomes(outcomes):
    # Deferred so that a failed call raises at its own position in the output
    # sequence, matching Executor.map().
    for ok, value in outcomes:
        if ok:
            yield value
        else:
            raise value


class SerialExecutor(Executor):
    """An :class:`~concurrent.futures.Executor` that runs each call immediately
    in the calling thread.

    No threads, processes, queues, or background state are created. A call
    submitted to this executor has finished by the time :meth:`submit` returns,
    so the :class:`~concurrent.futures.Future` it hands back is always complete
    and :meth:`Future.result` never blocks.

    Useful as a drop-in stand-in for a real executor when parallelism is
    unwanted (debugging, profiling, single-core runs, or tasks too small to be
    worth dispatching), and as a baseline when measuring the cost of
    parallelism.

    Examples
    --------

    >>> from westpa.work_managers.executors import SerialExecutor
    >>> with SerialExecutor() as executor:
    ...     future = executor.submit(pow, 2, 10)
    ...     future.done()
    ...     future.result()
    ...
    True
    1024

    Notes
    -----
    Differences from the pooled executors in :mod:`concurrent.futures`, all of
    them consequences of running inline:

    - Calls run in the thread that submits them, not a worker thread. A
      callable that inspects thread-local state or thread identity sees the
      caller's.
    - Nothing is ever queued, so a submitted call cannot be cancelled;
      :meth:`Future.cancel` always returns False, and the ``cancel_futures``
      argument to :meth:`shutdown` has nothing to act on.
    - Submitting from inside a running task is safe and cannot deadlock, unlike
      :class:`~concurrent.futures.ThreadPoolExecutor`, where a task that waits
      on another task can exhaust the pool.
    - :meth:`map` ignores ``timeout`` (results are always ready) and
      ``chunksize`` (there are no workers to distribute to).

    """

    def __init__(self):
        self._shutdown = False

    def submit(self, fn, /, *args, **kwargs):
        """Run ``fn(*args, **kwargs)`` and return a completed future.

        Parameters
        ----------
        fn : callable
            Callable to invoke.
        *args, **kwargs
            Arguments to pass to `fn`.

        Returns
        -------
        future : concurrent.futures.Future
            Future holding the return value of `fn`, or the exception it
            raised. The future is already done.

        Raises
        ------
        RuntimeError
            If the executor has been shut down. Exceptions raised by `fn`
            itself are stored on the future rather than propagated, as with the
            pooled executors.

        """
        if self._shutdown:
            raise RuntimeError('cannot schedule new futures after shutdown')

        future = Future()
        try:
            result = fn(*args, **kwargs)
        except BaseException as exc:
            future.set_exception(exc)
        else:
            future.set_result(result)
        return future

    def map(self, fn, *iterables, timeout=None, chunksize=1):
        """Apply `fn` to each set of arguments drawn from `iterables`.

        Every call is made before the returned iterator yields anything, as
        with :meth:`Executor.map`. A call that raises re-raises its exception
        at the corresponding position in the output.

        This path allocates no futures, so it carries no per-call overhead
        beyond calling `fn`.

        Returns
        -------
        results : iterator
            Iterator over the return values, in input order.

        """
        if self._shutdown:
            raise RuntimeError('cannot schedule new futures after shutdown')

        outcomes = []
        for args in zip(*iterables):
            try:
                outcomes.append((True, fn(*args)))
            except BaseException as exc:
                outcomes.append((False, exc))

        return _iter_outcomes(outcomes)

    def shutdown(self, wait=True, *, cancel_futures=False):
        """Reject further calls to :meth:`submit` and :meth:`map`.

        Nothing is pending when this is called, so `wait` and `cancel_futures`
        have no effect.

        """
        self._shutdown = True
