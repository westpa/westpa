import os, signal, tempfile, time

from work_managers import WMFuture
from work_managers.zeromq import ZMQWorkManager, ZMQWMMaster
from tsupport import *

import zmq

import nose.tools
from nose.tools import raises, nottest

class BaseTestZMQWorkManager:                
    @nose.tools.timed(2)
    def test_shutdown(self):
        work_manager = self.test_master
        work_manager.startup()
        #time.sleep(0.1)
        work_manager.shutdown()
        time.sleep(0.1)
        work_manager.remove_ipc_endpoints()
        assert not work_manager._dispatch_thread.is_alive(), 'dispatch thread still alive'
        assert not work_manager._ipc_endpoints, 'IPC endpoints not deleted'
        
    @nose.tools.timed(2)
    def test_submit(self):
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        task_socket = self.test_client_context.socket(zmq.PULL)
        task_socket.connect(task_endpoint)
        
        future = work_manager.submit(identity, 1)
        assert isinstance(future, WMFuture)
        task_object = task_socket.recv_pyobj()
        assert task_object[1] == identity
        assert tuple(task_object[2]) == (1,)
        task_socket.close()

    @nose.tools.timed(2)
    def test_multi_submit(self):
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        task_socket = self.test_client_context.socket(zmq.PULL)
        task_socket.connect(task_endpoint)
        
        futures = [work_manager.submit(identity, i) for i in xrange(5)]
        params = set()
        for i in xrange(5):
            (task_id, fn, args, kwargs) = task_socket.recv_pyobj()
            assert fn == identity
            params.add(args)
        assert params == set((i,) for i in xrange(5))
        task_socket.close()

class TestZMQWorkManagerIPC(BaseTestZMQWorkManager):
    def setUp(self):
        self.test_client_context = zmq.Context()
        task_endpoint = ZMQWorkManager.make_ipc_endpoint()
        result_endpoint = ZMQWorkManager.make_ipc_endpoint()
        ann_endpoint = ZMQWorkManager.make_ipc_endpoint()
        self.test_master = ZMQWMMaster(None, ann_endpoint, task_endpoint, result_endpoint)
        
    def tearDown(self):
        self.test_master.shutdown()
        self.test_master.remove_ipc_endpoints()
        self.test_client_context.term()
        del self.test_master, self.test_client_context

class TestZMQWorkManagerTCP(BaseTestZMQWorkManager):
    def setUp(self):
        self.test_client_context = zmq.Context()
        ann_endpoint = 'tcp://127.0.0.1:23811'
        task_endpoint = 'tcp://127.0.0.1:23812'
        result_endpoint = 'tcp://127.0.0.1:23813'
        self.test_master = ZMQWMMaster(None, ann_endpoint, task_endpoint, result_endpoint)

    def tearDown(self):
        self.test_master.shutdown()
        self.test_client_context.term()
        del self.test_master, self.test_client_context

            