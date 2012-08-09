from __future__ import division, print_function; __metaclass__ = type
import os, signal, tempfile, time, sys, multiprocessing, uuid
import cPickle as pickle

from work_managers import WMFuture
from work_managers.zeromq import Task, Result, ZMQBase, ZMQWorkManager, ZMQWMProcess, ZMQClient, recvall
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
            
            task = Task.from_zmq_frames(task_socket.recv_multipart(copy=False))
            assert task.server_id == work_manager.instance_id
            assert task.fn == identity
            assert task.args == (1,)
            
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
            for i in xrange(5):
                work_manager.submit(identity,i)

            params = set()
            for i in xrange(5):
                task = Task.from_zmq_frames(task_socket.recv_multipart(copy=False))
                assert task.server_id == work_manager.instance_id
                assert task.fn == identity
                params.add(task.args)
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
            task = Task.from_zmq_frames(task_socket.recv_multipart(copy=False))
            result = Result(task.server_id, task.task_id, value = 1)
            result_socket.send_multipart(result.to_zmq_frames())                        
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
        self.task_endpoint = 'tcp://127.0.0.1:23812'
        self.result_endpoint = 'tcp://127.0.0.1:23813'
                
        self.wmproc = ZMQWMProcess(self.task_endpoint, self.result_endpoint)
        self.wmproc.start()
        
        self.context = zmq.Context()
        
        self.task_socket = self.context.socket(zmq.REP)
        self.task_socket.bind(self.task_endpoint)
        
        self.result_socket = self.context.socket(zmq.PULL)
        self.result_socket.bind(self.result_endpoint)
        
        self.server_id = uuid.uuid4()
        self.node_id =  uuid.uuid4()
                

    def tearDown(self):
        self.context.destroy(linger=0)
        self.wmproc.terminate()
        del self.context, self.wmproc, self.server_id, self.node_id
    
    @nose.tools.timed(2)
    def test_task(self):
        task_id = uuid.uuid4()
        task = Task(self.server_id, task_id, identity, (1,), {})
        
        # admit request
        frames = self.task_socket.recv_multipart(copy=False)
        (node_id, client_id, message) = pickle.loads(frames[0].buffer.tobytes())
        assert node_id is None
        assert message == 'task'
        
        # reply to request for task
        self.task_socket.send_pyobj((self.node_id, client_id, 'task'), flags=zmq.SNDMORE)
        self.task_socket.send_multipart(task.to_zmq_frames())
        
        # receive result
        frames = self.result_socket.recv_multipart(copy=False)
        (node_id, client_id_2, message) = pickle.loads(frames[0].buffer.tobytes())
        
        
        assert node_id == self.node_id
        assert client_id_2 == client_id
        assert message == 'result'
        
        result = Result.from_zmq_frames(frames[1:])
        assert result.server_id == self.server_id
        assert result.task_id == task_id
        assert result.is_result
        assert result.value == 1
            
    @nose.tools.timed(2)
    def test_multiple_tasks(self):  
        
        result_values = set()
        for n in xrange(4):
            task_id = uuid.uuid4()
            task = Task(self.server_id, task_id, identity, (n,), {})
            
            # admit request
            frames = self.task_socket.recv_multipart(copy=False)
            (node_id, client_id, message) = pickle.loads(frames[0].buffer.tobytes())
            if n == 0:
                assert node_id is None
            else:
                assert node_id == self.node_id
            assert message == 'task'
            
            # reply to request for task
            self.task_socket.send_pyobj((self.node_id, client_id, 'task'), flags=zmq.SNDMORE)
            self.task_socket.send_multipart(task.to_zmq_frames())
            
            # receive result
            frames = self.result_socket.recv_multipart(copy=False)
            (node_id, client_id_2, message) = pickle.loads(frames[0].buffer.tobytes())
            
            
            assert node_id == self.node_id
            assert client_id_2 == client_id
            assert message == 'result'
            
            result = Result.from_zmq_frames(frames[1:])
            assert result.server_id == self.server_id
            assert result.task_id == task_id
            assert result.is_result
            result_values.add(result.value)
            
        assert result_values == set(xrange(4))

    @nose.tools.timed(2)
    def test_task_exception(self):
        task_id = uuid.uuid4()
        task = Task(self.server_id, task_id, will_fail, (), {})
        
        # admit request
        frames = self.task_socket.recv_multipart(copy=False)
        (node_id, client_id, message) = pickle.loads(frames[0].buffer.tobytes())
        assert node_id is None
        assert message == 'task'
        
        # reply to request
        self.task_socket.send_pyobj((self.node_id, client_id, 'task'), flags=zmq.SNDMORE)
        self.task_socket.send_multipart(task.to_zmq_frames())
        
        # receive result
        frames = self.result_socket.recv_multipart(copy=False)
        (node_id, client_id_2, message) = pickle.loads(frames[0].buffer.tobytes())
        
        
        assert node_id == self.node_id
        assert client_id_2 == client_id
        assert message == 'result'
        
        result = Result.from_zmq_frames(frames[1:])
        assert result.server_id == self.server_id
        assert result.task_id == task_id
        assert result.is_exception
        assert isinstance(result.exception, ExceptionForTest)
 
class TestZMQClient:
    def setUp(self):
        self.ann_endpoint = 'tcp://127.0.0.1:23811'
        self.task_endpoint = 'tcp://127.0.0.1:23812'
        self.result_endpoint = 'tcp://127.0.0.1:23813'
                
        self.server_id = uuid.uuid4()
        self.node_id =  uuid.uuid4()
        
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint, 4)
        
        self.context = zmq.Context()
                
    def tearDown(self):
        #self.test_client.shutdown_workers()
        self.context.destroy(linger=0)
        del self.context, self.server_id, self.node_id
    
    @nose.tools.timed(2)
    def test_spawn_workers(self):
        self.test_client.spawn_workers()
        time.sleep(0.1)
        
        try:
            for proc in self.test_client.workers:
                assert proc.is_alive()
        finally:        
            self.test_client.shutdown_workers()
        
    @nose.tools.timed(2)    
    def test_shutdown_workers(self):
        self.test_client.spawn_workers()
        time.sleep(0.1)
        self.test_client.shutdown_workers()
        time.sleep(0.1)
        assert not self.test_client.workers

    @timed(2)
    def test_startup(self):
        self.test_client.startup()
        
        try:
            assert self.test_client.workers
            assert self.test_client._taskfwd_thread is not None and self.test_client._taskfwd_thread.is_alive()
            assert self.test_client._rslfwd_thread is not None and self.test_client._rslfwd_thread.is_alive()
            assert self.test_client._monitor_thread is not None and self.test_client._monitor_thread.is_alive()
        finally:
            self.test_client.shutdown()
                 
    @timed(2)
    def test_shutdown(self):
        self.test_client.startup()
        time.sleep(0.1)
        self.test_client.shutdown()
        time.sleep(0.1)
        
        assert not self.test_client.workers
        assert self.test_client._monitor_thread is not None and not self.test_client._monitor_thread.is_alive()
        # cannot test task forwarder, as it deliberately blocks
        assert self.test_client._rslfwd_thread is not None and not self.test_client._rslfwd_thread.is_alive()
        
        
