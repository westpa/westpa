"""
An interface for using Ray as a work manager.
"""

import ray
import logging
from .core import WorkManager
from itertools import islice

log = logging.getLogger(__name__)


class RayFuture:
    """
    The way the different work managers were set up, they seem to generally use a pattern where a number of tasks
        are submitted, which return futures.
    The futures are looped over, and results returned.
    This duplicates some functionality in Ray, and I think at least for the Ray work manager, we don't need to use them,
        though we may well decide to continue using them anyways.

    For this demo, I'm going to try eliminating them for simplicity.
    However, I need to be sure the same functions are exposed so we don't have to change anything elsewhere.
    """

    def __init__(self, data):

        self.data = data

    def get_result(self, discard):
        """
        Return the result of a future.
        Although the discard argument is unused here, it's passed in some places, so this function needs to be able
            to accept it.
        """

        return self.data


class RayWorkManager(WorkManager):
    """
    A work manager using the ``ray`` module.
    """

    @classmethod
    def from_environ(cls, wmenv=None):
        """
        Returns an instance of the current work manager, set up with whatever parameters the environment stores.
        """

        return cls()

    def submit_as_completed(self, task_generator, queue_size=None):
        """
        Take in a generator of tasks.
        Create the tasks, send them to Ray, and return the results.
        """

        pending_ray_ids = set()

        for (fn, args, kwargs) in islice(task_generator, queue_size):

            ray_id = self.ray_submit(fn, args, kwargs)

            pending_ray_ids.add(ray_id)

        # Wait for all results to be completed, and spit them out
        results = ray.get(list(pending_ray_ids))
        result_futures = [RayFuture(result) for result in results]

        return result_futures

    def startup(self):
        """
        Initialize Ray.

        In the future, if running on a cluster, you could initialize Ray workers on a number of nodes.
            Then, this could just connect to that existing cluster.
            Example: https://github.com/jdrusso/ray_distributed
        A set of arguments to the work-manager could tell it whether to look for a cluster at a certain IP, or
            to just set up its own on the current node.
        """

        ray.init()

    def shutdown(self):
        """
        Cleanly close Ray.
        """

        self.shutdown_received = True
        ray.shutdown()

    def __init__(self, n_workers=None, shutdown_timeout=1):
        """
        Initialize the Ray work-manager.
        """

        super().__init__()

        self.n_workers = n_workers
        self.workers = None
        self.shutdown_received = False

    def ray_submit(self, fn, args=(), kwargs={}):
        """
        Accept a function, and farm it out to Ray workers. Does NOT wait for completion or return the result.

        This is *not* yet a drop-in replacement for submit(), which returns a Future.
        """

        @ray.remote
        def remote_function():
            """A generic wrapper to make a function remote and Ray-able."""
            return fn(*args, **kwargs)

        ray_id = remote_function.remote()

        return ray_id
