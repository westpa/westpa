from __future__ import division, print_function; __metaclass__ = type
import os, signal, tempfile, time, sys, multiprocessing, uuid, socket, json, re
import cPickle as pickle

from work_managers import WMFuture
from work_managers.zeromq import Task, Result, ZMQBase, ZMQWorkManager, ZMQWMProcess, ZMQClient, recvall, WorkerTerminated,\
    ZMQWMServer
from tsupport import *

import zmq

import nose.tools
from nose.tools import raises, nottest, timed

def sockdelay():
    '''Delay for slightly longer than the default auto-reconnect time for ZeroMQ (100 ms)'''
    time.sleep(0.11)
    
def randport():
    s = socket.socket()
    s.bind(('127.0.0.1',0))
    port = s.getsockname()[1]
    s.close()
    return port

def randipc():
    (fd, socket_path) = tempfile.mkstemp()
    os.close(fd)
    endpoint = 'ipc://{}'.format(socket_path)
    return endpoint
    

class BaseTestZMQWMServer:
                    
    @timed(2)
    def test_shutdown(self):
        work_manager = self.test_master
        work_manager.startup()
        work_manager.shutdown()
        sockdelay()
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
            sockdelay()
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
    
    @timed(30)
    def test_overkill_multi_submit(self):
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        task_socket = self.test_client_context.socket(zmq.PULL)
        task_socket.connect(task_endpoint)
        
        try:
            for i in xrange(10000):
                work_manager.submit(identity,i)

            params = set()
            for i in xrange(10000):
                print(i)
                task = Task.from_zmq_frames(task_socket.recv_multipart(copy=False))
                assert task.server_id == work_manager.instance_id
                assert task.fn == identity
                params.add(task.args)
            assert params == set((i,) for i in xrange(10000))
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

class TestZMQWMServerIPC(BaseTestZMQWMServer):
    def setUp(self):
        self.test_client_context = zmq.Context()
        task_endpoint = ZMQBase.make_ipc_endpoint()
        result_endpoint = ZMQBase.make_ipc_endpoint()
        ann_endpoint = ZMQBase.make_ipc_endpoint()
        self.test_master = ZMQWMServer(task_endpoint, result_endpoint, ann_endpoint,10)
        
    def tearDown(self):
        self.test_master.shutdown()
        self.test_master.remove_ipc_endpoints()
        del self.test_master, self.test_client_context

#@nose.SkipTest
class TestZMQWMServerTCP(BaseTestZMQWMServer):
    def setUp(self):
        self.test_client_context = zmq.Context()
        ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.test_master = ZMQWMServer(task_endpoint, result_endpoint, ann_endpoint, 10)

    def tearDown(self):
        self.test_master.shutdown()
        self.test_client_context.destroy(linger=0)
        del self.test_master, self.test_client_context
        
class TestZMQWMProcess:
    def setUp(self):
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
                
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
        self.wmproc.terminate()
        del self.context, self.wmproc, self.server_id, self.node_id
    
    @timed(2)
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
            
    @timed(2)
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

    @timed(2)
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
 
class BaseTestZMQClient:    
    @timed(2)
    def test_spawn_workers(self):
        #self.test_client.spawn_workers()
        self.test_client.startup()
        sockdelay()
        
        try:
            for (_pid,proc) in self.test_client.workers.items():
                assert proc.is_alive()
        finally:        
            self.test_client.shutdown()
        
    @timed(2)    
    def test_shutdown_workers(self):
        #self.test_client.spawn_workers()
        self.test_client.startup()
        sockdelay()
        self.test_client.shutdown()
        sockdelay()
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
    def test_shutdown_internal(self):
        self.test_client.startup()
        sockdelay()
        self.test_client.shutdown()
        sockdelay()
        
        assert not self.test_client.workers
        assert self.test_client._monitor_thread is not None and not self.test_client._monitor_thread.is_alive()
        assert self.test_client._taskfwd_thread is not None and not self.test_client._taskfwd_thread.is_alive()
        assert self.test_client._rslfwd_thread is not None and not self.test_client._rslfwd_thread.is_alive()
    
    @timed(2)    
    def test_shutdown_external(self):
        
        self.test_client.startup()
        sockdelay()
        
        ann_socket = self.context.socket(zmq.PUB)
        ann_socket.bind(self.ann_endpoint)
        # this delay is critical, as it allows the auto-reconnect logic to connect this socket to the
        # client; otherwise, an immediate send will result in discard of the shutdown message
        sockdelay() 
        ann_socket.send('shutdown')
        ann_socket.close()
        self.test_client._wait_for_shutdown()
        
        assert not self.test_client.workers
        assert self.test_client._monitor_thread is not None and not self.test_client._monitor_thread.is_alive()
        assert self.test_client._taskfwd_thread is not None and not self.test_client._taskfwd_thread.is_alive()
        assert self.test_client._rslfwd_thread is not None and not self.test_client._rslfwd_thread.is_alive()
        
    # the critical feature of this test is the @timed(2); if the wait hangs, then this test
    # times out
    @timed(2)
    def test_wait_shutdown(self):
        self.test_client.startup()
        sockdelay()
        self.test_client.shutdown()
        
        assert not self.test_client.workers
        assert self.test_client._monitor_thread is not None and not self.test_client._monitor_thread.is_alive()
        assert self.test_client._taskfwd_thread is not None and not self.test_client._taskfwd_thread.is_alive()
        assert self.test_client._rslfwd_thread is not None and not self.test_client._rslfwd_thread.is_alive()
                
    @timed(2)        
    def test_task_forward(self):
        outgoing_task_socket = self.context.socket(zmq.PUSH)
        outgoing_task_socket.setsockopt(zmq.HWM,1)
        outgoing_task_socket.bind(self.task_endpoint)
                
        self.test_client.startup(spawn_workers=False)
        
        fake_worker_task_socket = self.context.socket(zmq.REQ)
        fake_worker_task_socket.connect(self.test_client.worker_task_endpoint)
        
        try:
            # make task available as the server would
            task_id = uuid.uuid4()
            outgoing_task = Task(self.server_id, task_id, identity, (1,), {})
            outgoing_task_socket.send_multipart(outgoing_task.to_zmq_frames(),copy=False)
            
            # request task as worker would
            fake_worker_task_socket.send_pyobj((None, os.getpid(), 'task'))
            
            # receive task as worker would
            frames = fake_worker_task_socket.recv_multipart(copy=False)
            incoming_task = Task.from_zmq_frames(frames[1:])
            
            assert incoming_task.fn == identity
            assert incoming_task.args == (1,)
            assert incoming_task.kwargs == {}
            
        finally:
            outgoing_task_socket.close(linger=0)
            fake_worker_task_socket.close(linger=0)
            self.test_client.shutdown()

    @timed(2)            
    def test_result_forward(self):
        self.test_client.startup(spawn_workers=False)
        worker_result_socket = self.context.socket(zmq.PUSH)
        worker_result_socket.setsockopt(zmq.HWM,1)
        worker_result_socket.connect(self.test_client.worker_result_endpoint)
        
        master_result_socket = self.context.socket(zmq.PULL)
        master_result_socket.bind(self.test_client.upstream_result_endpoint)
        sockdelay()
        
        try:
            task_id = uuid.uuid4()
            outgoing_result = Result(self.server_id, task_id, value=1)
            worker_result_socket.send_pyobj((self.test_client.node_id, os.getpid(), 'result'), flags=zmq.SNDMORE)
            worker_result_socket.send_multipart(outgoing_result.to_zmq_frames(), copy=False)
            incoming_result = Result.from_zmq_frames(master_result_socket.recv_multipart(copy=False))
            
            print('outgoing = {}'.format(outgoing_result))
            print('incoming = {}'.format(incoming_result))

            assert incoming_result.server_id == outgoing_result.server_id
            assert incoming_result.task_id == outgoing_result.task_id
            assert incoming_result.value == outgoing_result.value
            
        finally:
            master_result_socket.close(linger=0)
            worker_result_socket.close(linger=0)
            self.test_client.shutdown()
    
    def test_dispatch(self):
        master_task_socket = self.context.socket(zmq.PUSH)
        master_task_socket.bind(self.task_endpoint)
        
        master_result_socket = self.context.socket(zmq.PULL)
        master_result_socket.bind(self.result_endpoint)
        
        self.test_client.startup()
        
        try:
            task_id = uuid.uuid4()
            task = Task(self.server_id, task_id, identity, (1,), {})
             
            master_task_socket.send_multipart(task.to_zmq_frames())
            
            result = Result.from_zmq_frames(master_result_socket.recv_multipart(copy=False))
            assert result.value == 1
            
        finally:
            master_result_socket.close()
            master_task_socket.close()
            self.test_client.shutdown()
            
    @timed(2)
    def test_dispatch_multiple(self):
        master_task_socket = self.context.socket(zmq.PUSH)
        master_task_socket.bind(self.task_endpoint)
        
        master_result_socket = self.context.socket(zmq.PULL)
        master_result_socket.bind(self.result_endpoint)
        
        self.test_client.startup()
        
        
        try:
            tasks = [Task(self.server_id, task_id=uuid.uuid4(), fn=identity, args=(n,), kwargs={})
                     for n in xrange(4)]
            
            for task in tasks:
                master_task_socket.send_multipart(task.to_zmq_frames(),copy=False)
                
            result_values = set()
            for _n in xrange(4):
                result = Result.from_zmq_frames(master_result_socket.recv_multipart(copy=False))
                result_values.add(result.value)
                
            assert result_values == set(xrange(4))
            
        finally:
            master_result_socket.close()
            master_task_socket.close()
            self.test_client.shutdown()
                    
    @timed(2)
    def test_missing_server_shutdown(self):
        announce_socket = self.context.socket(zmq.PUB)
        announce_socket.bind(self.ann_endpoint)
        
        try:
            self.test_client.server_heartbeat_interval = 0.1
            self.test_client.startup()
            sockdelay()
            announce_socket.send('ping')
            sockdelay()
        finally:
            announce_socket.close(linger=0)
        
    @timed(5)
    def test_autokill_hung_worker(self):        
        master_task_socket = self.context.socket(zmq.PUSH)
        master_task_socket.bind(self.task_endpoint)
        
        master_result_socket = self.context.socket(zmq.PULL)
        master_result_socket.bind(self.result_endpoint)
        
        self.test_client.worker_task_timeout = 0.2
        self.test_client.hangcheck = 0.05
        self.test_client.startup()
        
        try:
            # dispatch a pathological task
            task = Task(self.server_id, uuid.uuid4(), will_busyhang_uninterruptible, (), {})
            master_task_socket.send_multipart(task.to_zmq_frames(), copy=False)
            
            result = Result.from_zmq_frames(master_result_socket.recv_multipart(copy=False,))
            assert isinstance(result.exception, WorkerTerminated)
            
            self.test_client.shutdown()
            
        finally:
            master_result_socket.close()
            master_task_socket.close()
            
    

#@nose.SkipTest            
class TestZMQClientTCP(BaseTestZMQClient):
    def setUp(self):
        self.ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
                
        self.server_id = uuid.uuid4()
        self.node_id =  uuid.uuid4()
        
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint, 2)
        
        self.context = zmq.Context()
                
    def tearDown(self):
        #self.test_client.shutdown_workers()
        if self.test_client._monitor_thread is not None and self.test_client._monitor_thread.is_alive():
            self.test_client.shutdown()
        del self.context, self.server_id, self.node_id

class TestZMQClientIPC(BaseTestZMQClient):
    def setUp(self):
        self.ann_endpoint = ZMQBase.make_ipc_endpoint()
        self.task_endpoint = ZMQBase.make_ipc_endpoint()
        self.result_endpoint = ZMQBase.make_ipc_endpoint()
                
        self.server_id = uuid.uuid4()
        self.node_id =  uuid.uuid4()
        
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint, 2)
        
        self.context = zmq.Context()
                
    def tearDown(self):
        del self.context, self.server_id, self.node_id
        
class BaseTestCoordinated(CommonParallelTests,CommonWorkManagerTests):
    
    @timed(2)
    def test_sigint_shutdown(self):
        ann_socket = self.context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')        
        work_manager = self.work_manager
        work_manager.install_sigint_handler()
        
        ann_socket.connect(work_manager.master_announce_endpoint)
        sockdelay()
    
        try:
            os.kill(os.getpid(), signal.SIGINT)
        except KeyboardInterrupt:
            sockdelay()
            announcements = [ann_socket.recv()]
            announcements.extend(recvall(ann_socket))
            assert 'shutdown' in announcements            
            assert not work_manager._dispatch_thread.is_alive(), 'dispatch thread still alive'
            assert not work_manager._receive_thread.is_alive(), 'receive thread still alive'
            assert not work_manager._announce_thread.is_alive(), 'announcement thread still alive'
        finally:
            ann_socket.close(linger=0)
   
                
    
class TestCoordinatedIPC(BaseTestCoordinated):
    def setUp(self):
        self.nprocs = 4        
        
        self.ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        
        self.test_master = ZMQWorkManager(self.nprocs, self.task_endpoint, self.result_endpoint, self.ann_endpoint)
        self.test_master.startup()
        
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint, self.nprocs)
        self.test_client.startup()
        
        self.work_manager = self.test_master
        
        self.context = zmq.Context()
    
    def tearDown(self):
        self.test_client.shutdown()
        self.test_master.shutdown()
        
class TestCoordinatedTCP(BaseTestCoordinated):
    def setUp(self):
        self.nprocs = 4
        
        self.ann_endpoint = ZMQBase.make_ipc_endpoint()
        self.task_endpoint = ZMQBase.make_ipc_endpoint()
        self.result_endpoint = ZMQBase.make_ipc_endpoint()
        
        self.test_master = ZMQWorkManager(self.nprocs, self.task_endpoint, self.result_endpoint, self.ann_endpoint)
        self.test_master.startup()
        
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint, self.nprocs)
        self.test_client.startup()
        
        self.work_manager = self.test_master        

        self.context = zmq.Context()
    
    def tearDown(self):
        self.test_client.shutdown()
        self.test_master.shutdown()
        
class TestZMQWorkManager:
    sanitize_vars = ('WM_WORK_MANAGER', 'WM_N_WORKERS',
                     'WM_ZMQ_SERVER_INFO', 'WM_ZMQ_COMM_MODE',
                     'WM_ZMQ_TASK_ENDPOINT', 'WM_ZMQ_RESULT_ENDPOINT', 'WM_ZMQ_ANNOUNCE_ENDPOINT')
    
    def setUp(self):
        for varname in self.sanitize_vars:
            assert varname not in os.environ
    
    def tearDown(self):
        for varname in self.sanitize_vars:
            os.environ.pop(varname, None)
                    
    def test_auto_local(self):
        with ZMQWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
            
    def test_server_info_ipc(self):
        with ZMQWorkManager() as work_manager:
            server_info_filename = work_manager.server_info_filename
            assert os.path.exists(server_info_filename)
            assert os.stat(server_info_filename).st_mode & 0777 == 0600
            server_info = json.load(open(server_info_filename, 'rt'))
            assert re.sub(r'\*', socket.gethostname(), work_manager.master_task_endpoint) == server_info['task_endpoint']
            assert re.sub(r'\*', socket.gethostname(), work_manager.master_result_endpoint) == server_info['result_endpoint']
            assert re.sub(r'\*', socket.gethostname(), work_manager.master_announce_endpoint) == server_info['announce_endpoint']
        assert not os.path.exists(server_info_filename)
        
    def test_environ_empty(self):
        with ZMQWorkManager.from_environ() as work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
    
    def test_environ_mode_tcp(self):
        os.environ['WM_ZMQ_COMM_MODE'] = 'tcp'
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.master_task_endpoint.startswith('tcp://*')
            assert work_manager.master_result_endpoint.startswith('tcp://*')
            assert work_manager.master_announce_endpoint.startswith('tcp://*')
            future = work_manager.submit(will_succeed)
            future.get_result()

    def test_environ_ipc_endpoints(self):
        task_endpoint = randipc()
        result_endpoint = randipc()
        announce_endpoint = randipc()
        
        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.master_task_endpoint == task_endpoint
            assert work_manager.master_result_endpoint == result_endpoint
            assert work_manager.master_announce_endpoint == announce_endpoint
            future = work_manager.submit(will_succeed)
            future.get_result()

    def test_environ_tcp_endpoints(self):
        # note that this tests not only that the work manager honor our environment settings, but that
        # the hostname-to-ip mapping succeeded
        task_endpoint = 'tcp://localhost:{}'.format(randport())
        result_endpoint = 'tcp://localhost:{}'.format(randport())
        announce_endpoint = 'tcp://localhost:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.master_task_endpoint == re.sub('localhost','127.0.0.1', task_endpoint)
            assert work_manager.master_result_endpoint == re.sub('localhost','127.0.0.1', result_endpoint)
            assert work_manager.master_announce_endpoint == re.sub('localhost','127.0.0.1', announce_endpoint)
            future = work_manager.submit(will_succeed)
            future.get_result()

    def test_environ_nworkers(self):
        os.environ['WM_N_WORKERS'] = str(2)
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.internal_client.n_workers == 2
            future = work_manager.submit(will_succeed)
            future.get_result()

    def test_environ_noworkers(self):
        os.environ['WM_N_WORKERS'] = str(0)
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.internal_client is None

    @timed(2)
    @raises(ValueError)
    def test_client_from_bad_environ(self):
        os.environ['WM_N_WORKERS'] = str(0)
        task_endpoint = 'tcp://*:{}'.format(randport())
        result_endpoint = 'tcp://*:{}'.format(randport())
        announce_endpoint = 'tcp://*:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        
        with ZMQWorkManager() as work_manager:
            os.environ['WM_N_WORKERS'] = str(2)
            test_client = ZMQClient.from_environ()
            test_client.startup()            
            try:
                future = work_manager.submit(will_succeed)
                future.get_result()                
            finally:
                test_client.shutdown()
                
    @timed(2)
    def test_client_from_environ(self):
        os.environ['WM_N_WORKERS'] = str(0)
        task_endpoint = 'tcp://localhost:{}'.format(randport())
        result_endpoint = 'tcp://localhost:{}'.format(randport())
        announce_endpoint = 'tcp://localhost:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        
        with ZMQWorkManager() as work_manager:
            os.environ['WM_N_WORKERS'] = str(2)
            test_client = ZMQClient.from_environ()
            test_client.startup()            
            try:
                future = work_manager.submit(will_succeed)
                future.get_result()                
            finally:
                test_client.shutdown()
                
    @nose.tools.timed(2)                    
    def test_worker_ids(self):
        work_manager = ZMQWorkManager()
        with work_manager:
            futures = work_manager.submit_many([(get_process_index, (), {})] * work_manager.n_workers)
            work_manager.wait_all(futures)
            results = set(future.get_result() for future in futures)
            assert results == set(str(n) for n in xrange(work_manager.n_workers))
        