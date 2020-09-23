import itertools
import time
import unittest

from westpa.work_managers.zeromq import ZMQWorkManager, ZMQWorker
from westpa.work_managers.zeromq.node import ZMQNode

from ..tsupport import CommonWorkManagerTests
from ..tsupport import identity, random_int

from .zmq_tsupport import SETUP_WAIT, TEARDOWN_WAIT, SHUTDOWN_WAIT, BEACON_PERIOD, BEACON_WAIT
from .zmq_tsupport import ZMQTestBase


class TestZMQNodeExternal(ZMQTestBase, CommonWorkManagerTests, unittest.TestCase):
    n_workers = 2

    '''Tests for the core task dispersal/retrieval and shutdown operations
    (the parts of the WM that do not require ZMQWorker).'''

    def setUp(self):
        super().setUp()

        self.test_wm = ZMQWorkManager(n_local_workers=0)
        upstream_ann_endpoint = self.test_core.make_internal_endpoint()
        upstream_rr_endpoint = self.test_core.make_internal_endpoint()
        downstream_ann_endpoint = self.test_core.make_internal_endpoint()
        downstream_rr_endpoint = self.test_core.make_internal_endpoint()
        self.test_wm.downstream_rr_endpoint = upstream_rr_endpoint
        self.test_wm.downstream_ann_endpoint = upstream_ann_endpoint

        self.test_node = ZMQNode(upstream_rr_endpoint, upstream_ann_endpoint, 0)
        self.test_node.downstream_ann_endpoint = downstream_ann_endpoint
        self.test_node.downstream_rr_endpoint = downstream_rr_endpoint

        self.test_workers = [ZMQWorker(downstream_rr_endpoint, downstream_ann_endpoint) for n in range(self.n_workers)]

        # Set operation parameters
        for core_object in itertools.chain([self.test_wm, self.test_node, self.test_core], self.test_workers):
            core_object.validation_fail_action = 'raise'
            core_object.master_beacon_period = BEACON_PERIOD
            core_object.task_beacon_period = BEACON_PERIOD

        for worker in self.test_workers:
            worker.master_beacon_period = BEACON_WAIT
            worker.shutdown_timeout = 0.5
            worker.startup()

        self.test_node.startup()
        self.test_wm.startup()

        # self.test_wm.startup()
        # self.test_node.startup()
        # self.test_worker.startup()

        self.test_core.master_id = self.test_wm.master_id
        self.test_node.master_id = self.test_wm.master_id

        self.work_manager = self.test_wm

        # time.sleep(1.0)
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        self.test_wm.signal_shutdown()
        time.sleep(TEARDOWN_WAIT)

        self.test_wm.comm_thread.join()

        self.test_node.signal_shutdown()
        self.test_node.comm_thread.join()

        for worker in self.test_workers:
            worker.signal_shutdown()
            worker.comm_thread.join()

        self.test_node.remove_ipc_endpoints()
        self.test_wm.remove_ipc_endpoints()

        super().tearDown()

    def test_shutdown(self):
        time.sleep(SHUTDOWN_WAIT)
        self.test_wm.signal_shutdown()
        self.test_node.join()
        assert not self.test_node.comm_thread.is_alive()

    def test_task(self):
        r = random_int()
        future = self.test_wm.submit(identity, (r,), {})
        assert future.get_result() == r


class TestZMQNodeInternal(ZMQTestBase, CommonWorkManagerTests, unittest.TestCase):
    n_workers = 2

    '''Tests for the core task dispersal/retrieval and shutdown operations
    (the parts of the WM that do not require ZMQWorker).'''

    def setUp(self):
        super().setUp()

        self.test_wm = ZMQWorkManager(n_local_workers=0)
        upstream_ann_endpoint = self.test_core.make_internal_endpoint()
        upstream_rr_endpoint = self.test_core.make_internal_endpoint()
        downstream_ann_endpoint = self.test_core.make_internal_endpoint()
        downstream_rr_endpoint = self.test_core.make_internal_endpoint()
        self.test_wm.downstream_rr_endpoint = upstream_rr_endpoint
        self.test_wm.downstream_ann_endpoint = upstream_ann_endpoint

        self.test_node = ZMQNode(upstream_rr_endpoint, upstream_ann_endpoint, self.n_workers)
        self.test_node.downstream_ann_endpoint = downstream_ann_endpoint
        self.test_node.downstream_rr_endpoint = downstream_rr_endpoint

        # Set operation parameters
        for core_object in itertools.chain([self.test_wm, self.test_node, self.test_core]):
            core_object.validation_fail_action = 'raise'
            core_object.master_beacon_period = BEACON_PERIOD
            core_object.task_beacon_period = BEACON_PERIOD

        for worker in self.test_node.local_workers:
            worker.master_beacon_period = BEACON_WAIT
            worker.shutdown_timeout = 0.5

        self.test_node.startup()
        self.test_wm.startup()

        # self.test_wm.startup()
        # self.test_node.startup()
        # self.test_worker.startup()

        self.test_core.master_id = self.test_wm.master_id
        self.test_node.master_id = self.test_wm.master_id

        self.work_manager = self.test_wm

        # time.sleep(1.0)
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        self.test_wm.signal_shutdown()
        time.sleep(TEARDOWN_WAIT)

        self.test_wm.comm_thread.join()

        self.test_node.signal_shutdown()
        self.test_node.comm_thread.join()

        self.test_node.remove_ipc_endpoints()
        self.test_wm.remove_ipc_endpoints()

        super().tearDown()
