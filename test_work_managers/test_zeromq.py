from __future__ import division, print_function; __metaclass__ = type
import os, signal, tempfile, time, sys, multiprocessing, uuid

from work_managers import WMFuture
from work_managers.zeromq import ZMQBase, ZMQWorkManager, ZMQWMProcess, recvall
from tsupport import *

import zmq

import nose.tools
from nose.tools import raises, nottest, timed

class BaseTestZMQWMServer:
                    
    @timed(2)
    def test_shutdown(self):
        work_manager = self.test_master
        work_manager.startup()
        work_manager.shutdown()
        time.sleep(0.1)
        work_manager.remove_ipc_endpoints()
        assert not work_manager._dispatch_thread.is_alive(), 'dispatch thread still alive'
        assert not work_manager._receive_thread.is_alive(), 'receive thread still alive'
        assert not work_manager._announce_thread.is_alive(), 'announcement thread still alive'
        assert not work_manager._ipc_endpoints, 'IPC endpoints not deleted'
        
    @timed(2)
    def test_shutdown_announced(self):
        work_manager = self.test_master
        ann_socket = self.test_client_context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')        
        
        try:
            work_manager.startup()
            ann_socket.connect(work_manager.master_announce_endpoint)
            time.sleep(0.1)
            work_manager.shutdown()
            
            #announcements = recvall(ann_socket)
            ann = ann_socket.recv()
            #assert 'shutdown' in announcements
            assert ann == 'shutdown'
        finally:
            ann_socket.close(linger=0)
        
    @timed(2)
    def test_submit(self):
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        task_socket = self.test_client_context.socket(zmq.PULL)
        task_socket.connect(task_endpoint)
        
        try:
            future = work_manager.submit(identity, 1)
            assert isinstance(future, WMFuture)
            task_object = task_socket.recv_pyobj()
            assert task_object[0] == work_manager.instance_id
            assert task_object[2] == identity
            assert tuple(task_object[3]) == (1,)
        finally:
            task_socket.close(linger=0)

    @timed(2)
    def test_multi_submit(self):
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        task_socket = self.test_client_context.socket(zmq.PULL)
        task_socket.connect(task_endpoint)
        
        try:
            futures = [work_manager.submit(identity, i) for i in xrange(5)]
            params = set()
            for i in xrange(5):
                (instance_id, task_id, fn, args, kwargs) = task_socket.recv_pyobj()
                assert fn == identity
                params.add(args)
            assert params == set((i,) for i in xrange(5))
        finally:
            task_socket.close(linger=0)
    
    @timed(2)
    def test_receive_result(self):
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        result_endpoint = work_manager.master_result_endpoint
        task_socket = self.test_client_context.socket(zmq.PULL)
        result_socket = self.test_client_context.socket(zmq.PUSH)
        task_socket.connect(task_endpoint)
        result_socket.connect(result_endpoint)
        
        try:
            future = work_manager.submit(will_succeed)
            
            (instance_id, task_id, fn, args, kwargs) = task_socket.recv_pyobj()
            result_tuple = (instance_id, task_id, 'result', 1)
            result_socket.send_pyobj(result_tuple)
            
            assert future.get_result() == 1
        finally:        
            task_socket.close(linger=0)
            result_socket.close(linger=0)
        
    @timed(2)
    def test_server_heartbeat(self):
        work_manager = self.test_master
        work_manager.server_heartbeat_interval = 0.1
        ann_socket = self.test_client_context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')        
        
        try:
            work_manager.startup()
            ann_socket.connect(work_manager.master_announce_endpoint)
            time.sleep(0.2)
            
            ann = ann_socket.recv()
            assert ann == 'ping'
            work_manager.shutdown()
            
        finally:
            ann_socket.close(linger=0)

    @nose.tools.timed(2)
    def test_sigint_shutdown(self):
        ann_socket = self.test_client_context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')        
        work_manager = self.test_master
        work_manager.install_sigint_handler()
        
        work_manager.startup()
        ann_socket.connect(work_manager.master_announce_endpoint)
        time.sleep(0.1)
    
        try:
            os.kill(os.getpid(), signal.SIGINT)
        except KeyboardInterrupt:
            time.sleep(0.1)
            announcements = [ann_socket.recv()]
            announcements.extend(recvall(ann_socket))
            assert 'shutdown' in announcements            
            assert not work_manager._dispatch_thread.is_alive(), 'dispatch thread still alive'
            assert not work_manager._receive_thread.is_alive(), 'receive thread still alive'
            assert not work_manager._announce_thread.is_alive(), 'announcement thread still alive'
        finally:
            ann_socket.close(linger=0)

class TestZMQWMServerIPC(BaseTestZMQWMServer):
    def setUp(self):
        self.test_client_context = zmq.Context()
        task_endpoint = ZMQBase.make_ipc_endpoint()
        result_endpoint = ZMQBase.make_ipc_endpoint()
        ann_endpoint = ZMQBase.make_ipc_endpoint()
        self.test_master = ZMQWorkManager(task_endpoint, result_endpoint, ann_endpoint)
        
    def tearDown(self):
        self.test_master.shutdown()
        self.test_master.remove_ipc_endpoints()
        self.test_client_context.destroy(linger=0)
        del self.test_master, self.test_client_context

#@nose.SkipTest
class TestZMQWMServerTCP(BaseTestZMQWMServer):
    def setUp(self):
        self.test_client_context = zmq.Context()
        ann_endpoint = 'tcp://127.0.0.1:23811'
        task_endpoint = 'tcp://127.0.0.1:23812'
        result_endpoint = 'tcp://127.0.0.1:23813'
        self.test_master = ZMQWorkManager(task_endpoint, result_endpoint, ann_endpoint)

    def tearDown(self):
        self.test_master.shutdown()
        self.test_client_context.destroy(linger=0)
        del self.test_master, self.test_client_context
        
class TestZMQWMProcess:
    def setUp(self):
        self.context = zmq.Context()
        self.task_endpoint = 'tcp://127.0.0.1:23812'
        self.result_endpoint = 'tcp://127.0.0.1:23813'
        
        self.task_socket = self.context.socket(zmq.PUSH)
        self.result_socket = self.context.socket(zmq.PULL)
        
        self.task_socket.bind(self.task_endpoint)
        self.result_socket.bind(self.result_endpoint)
        
        self.wmproc = ZMQWMProcess(self.task_endpoint,self.result_endpoint)
        self.wmproc.start()

    def tearDown(self):
        self.context.destroy(linger=0)
        self.wmproc.terminate()
        del self.context, self.wmproc
    
    @nose.tools.timed(2)
    def test_task(self):        
        server_id = uuid.uuid4()
        task_id = uuid.uuid4()
                
        task_tuple = (server_id, task_id, identity, (1,), {})
        self.task_socket.send_pyobj(task_tuple)
        result_tuple = self.result_socket.recv_pyobj()
        assert result_tuple == (server_id, task_id, 'result', 1)

    @nose.tools.timed(2)
    def test_failed_task(self):
        server_id = uuid.uuid4()
        task_id = uuid.uuid4()
                
        task_tuple = (server_id, task_id, will_fail, (), {})
        self.task_socket.send_pyobj(task_tuple)
        result_tuple = self.result_socket.recv_pyobj()
        assert result_tuple[:3] == (server_id, task_id, 'exception')
        assert isinstance(result_tuple[3][0], ExceptionForTest)
        
    @nose.tools.timed(2)
    def test_multiple_tasks(self):        
        server_id = uuid.uuid4()
        
        for _n in xrange(4):
            task_id = uuid.uuid4()
            task_tuple = (server_id, task_id, identity, (1,), {})
            self.task_socket.send_pyobj(task_tuple)
            result_tuple = self.result_socket.recv_pyobj()
            assert result_tuple == (server_id, task_id, 'result', 1)
                
#    @nose.tools.timed(2)
#    @nose.tools.raises(ValueError)
#    def test_server_mismatch(self):
#        server_id = uuid.uuid4()
#        bad_server_id = uuid.uuid4()
#        task_id = uuid.uuid4()
#                
#        task_tuple = (server_id, task_id, identity, (1,), {})
#        self.task_socket.send_pyobj(task_tuple)
#        result_tuple = self.result_socket.recv_pyobj()
#        assert result_tuple == (server_id, task_id, 'result', 1)
#                
#        task_tuple = (bad_server_id, task_id, identity, (1,), {})
#        self.task_socket.send_pyobj(task_tuple)
        
        
