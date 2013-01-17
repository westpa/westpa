from __future__ import division, print_function; __metaclass__ = type
import os, signal, tempfile, time, sys, multiprocessing, uuid, socket, json, re
import cPickle as pickle

from work_managers import WMFuture
from work_managers import zeromq as zwm
from work_managers.zeromq import Task, ZMQBase, ZMQWorkManager, ZMQWMProcess, ZMQClient, recvall, WorkerTerminated,\
    ZMQWMServer
from tsupport import *

import zmq

import nose.tools
from nose.tools import raises, nottest, timed, assert_raises #@UnresolvedImport
from nose.plugins.skip import SkipTest

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

class TestZMQWMServer:
    def setUp(self):
        self.test_client_context = zmq.Context()
        ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.test_master = ZMQWMServer(task_endpoint, result_endpoint, ann_endpoint, 10)

    def tearDown(self):
        self.test_master.shutdown()
        self.test_master._close_signal_sockets()
        self.test_client_context.destroy(linger=0)
        del self.test_master, self.test_client_context
    
    def send_task_request(self, socket):
        socket.send_pyobj((zwm.MSG_TASK_REQUEST, None, None, None))
        
    def receive_task(self, socket):
        frames = socket.recv_multipart()
        #(tag, server_id, client_id, task_id) = zwm.unpickle_frame(frames[0])
        #(fn, args, kwargs) = zwm.unpickle_frame(frames[1])
        return zwm.unpickle_frame(frames[1])

    def test_shutdown(self):
        '''Server: shutdown works properly'''
        work_manager = self.test_master
        work_manager.startup()
        work_manager.shutdown()
        sockdelay()
        work_manager.remove_ipc_endpoints()
        assert not work_manager._dispatch_thread.is_alive(), 'dispatch thread still alive'
        assert not work_manager._receive_thread.is_alive(), 'receive thread still alive'
        assert not work_manager._announce_thread.is_alive(), 'announcement thread still alive'
        assert not work_manager._ipc_endpoints, 'IPC endpoints not deleted'
        
    def test_shutdown_announced(self):
        '''Server: announces shutdown to clients'''
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
            assert ann == zwm.MSG_SHUTDOWN
        finally:
            ann_socket.close(linger=0)
            
    def test_notasks(self):
        '''Server: Request to a server with no enqueued tasks yields a "no tasks" message'''
        work_manager = self.test_master
        work_manager.task_queue_wait = 0.001
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        task_socket = self.test_client_context.socket(zmq.REQ)
        task_socket.connect(task_endpoint)
        
        try:
            self.send_task_request(task_socket)
            frames = task_socket.recv_multipart()
            (tag, server_id, _, _) = zwm.unpickle_frame(frames[0])
            assert tag == zwm.MSG_TASK_UNAVAILABLE
            assert server_id == work_manager.instance_id

        finally:
            task_socket.close(linger=0)        
        
    def test_submit(self):
        '''Server: Calling submit yields properly-formatted message'''
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        task_socket = self.test_client_context.socket(zmq.REQ)
        task_socket.connect(task_endpoint)
        
        try:
            future = work_manager.submit(identity, (1,), {})
            assert isinstance(future, WMFuture)
            self.send_task_request(task_socket)
            frames = task_socket.recv_multipart()
            (tag, server_id, _, task_id) = zwm.unpickle_frame(frames[0])
            (fn, args, kwargs) = zwm.unpickle_frame(frames[1])
            
            assert task_id is not None
            assert tag == zwm.MSG_TASK_AVAILABLE
            assert server_id == work_manager.instance_id
            assert fn == identity
            assert args == (1,)
            assert kwargs == {}
            
        finally:
            task_socket.close(linger=0)
 
    def test_receive_result(self):
        '''Server: Receives results successfully'''
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        result_endpoint = work_manager.master_result_endpoint
        task_socket = self.test_client_context.socket(zmq.REQ)
        result_socket = self.test_client_context.socket(zmq.REQ)
        task_socket.connect(task_endpoint)
        result_socket.connect(result_endpoint)
        
        try:
            future = work_manager.submit(identity)
            self.send_task_request(task_socket)
            frames = task_socket.recv_multipart()
            (_, server_id, client_id, task_id) = zwm.unpickle_frame(frames[0])
            result_socket.send_pyobj((zwm.MSG_RESULT_SUBMISSION, server_id, client_id, task_id), flags=zmq.SNDMORE)
            result_socket.send_pyobj((zwm.RESULT_TYPE_RETVAL, 1))            
            assert future.get_result() == 1
        finally:        
            task_socket.close(linger=0)
            result_socket.close(linger=0)
            
    @raises(ValueError)
    def test_receive_exception(self):
        '''Server: Receives exceptions successfully'''
        work_manager = self.test_master
        work_manager.startup()
        task_endpoint = work_manager.master_task_endpoint
        result_endpoint = work_manager.master_result_endpoint
        task_socket = self.test_client_context.socket(zmq.REQ)
        result_socket = self.test_client_context.socket(zmq.REQ)
        task_socket.connect(task_endpoint)
        result_socket.connect(result_endpoint)
        
        try:
            future = work_manager.submit(identity)
            self.send_task_request(task_socket)
            frames = task_socket.recv_multipart()
            (_, server_id, client_id, task_id) = zwm.unpickle_frame(frames[0])
            result_socket.send_pyobj((zwm.MSG_RESULT_SUBMISSION, server_id, client_id, task_id), flags=zmq.SNDMORE)
            result_socket.send_pyobj((zwm.RESULT_TYPE_EXCEPTION, (ValueError(42), 'traceback')))
            future.get_result()
            
        finally:        
            task_socket.close(linger=0)
            result_socket.close(linger=0)

    def test_server_heartbeat(self):
        '''Server: emits heartbeats'''
        work_manager = self.test_master
        work_manager.server_heartbeat_interval = 0.1
        ann_socket = self.test_client_context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')        
        
        try:
            work_manager.startup()
            ann_socket.connect(work_manager.master_announce_endpoint)
            time.sleep(0.2)
            
            ann = ann_socket.recv()
            assert ann == zwm.MSG_PING
            work_manager.shutdown()
            
        finally:
            ann_socket.close(linger=0)
        
class TestZMQWMProcess:
    def setUp(self):
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
                
        self.wmproc = ZMQWMProcess(self.task_endpoint, self.result_endpoint)
        self.wmproc.start()
        
        self.context = zmq.Context()
        
        self.task_socket = self.context.socket(zmq.REP)
        self.task_socket.bind(self.task_endpoint)
        
        self.result_socket = self.context.socket(zmq.REP)
        self.result_socket.bind(self.result_endpoint)
        
        self.server_id = uuid.uuid4()
        self.node_id =  uuid.uuid4()

    def tearDown(self):
        self.wmproc.terminate()
        self.wmproc._close_signal_sockets()
        del self.context, self.wmproc, self.server_id, self.node_id
    
    def reply_with_task(self, fn, args, kwargs, socket=None):
        socket = socket or self.task_socket
        task_id = uuid.uuid4()
        frames = socket.recv_multipart(copy=False)
        (node_tag, node_id, client_id, client_pid) = zwm.unpickle_frame(frames[0])
        (server_tag, server_id, client_id, payload) = zwm.unpickle_frame(frames[1])
        
        socket.send_pyobj((zwm.MSG_TASK_AVAILABLE, self.node_id, client_id, client_pid), flags=zmq.SNDMORE)
        socket.send_pyobj((zwm.MSG_TASK_AVAILABLE, self.server_id, client_id, task_id), flags=zmq.SNDMORE)
        socket.send_pyobj((fn, args, kwargs))
        
    def receive_result(self, socket=None):
        socket = socket or self.result_socket
        frames = socket.recv_multipart()
        socket.send_pyobj((zwm.MSG_ACK, None, None, None))
        
        (node_tag, node_id, client_id, client_pid, task_id) = zwm.unpickle_frame(frames[0])
        (server_tag, server_id, client_id, task_id) = zwm.unpickle_frame(frames[1])
        (result_type, payload) = zwm.unpickle_frame(frames[2])
        return (result_type, payload)
        
    def test_task(self):
        '''Process: requests and processes single successful task'''
        
        self.reply_with_task(identity, (1,), {})
        result_type, payload = self.receive_result()
        
        assert result_type == zwm.RESULT_TYPE_RETVAL
        assert payload == 1
        
    def test_exception(self):
        '''Process: requests and processes single failing task'''
        
        self.reply_with_task(will_fail, (), {})
        result_type, payload = self.receive_result()
        
        assert result_type == zwm.RESULT_TYPE_EXCEPTION
        (exception, traceback) = payload
        assert isinstance(exception, ExceptionForTest)
        
    def test_multiple_tasks(self):
        '''Process: requests and processes multiple tasks'''
        results = set()
        for n in xrange(20):
            self.reply_with_task(identity, (n,), {})
            result_type, payload = self.receive_result()
            results.add(payload)
            
        assert results == set(xrange(20))
        
class TestZMQClient:    
    def setUp(self):
        self.ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
                
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint, 2)
        self.node_id = self.test_client.instance_id
        self.server_id = uuid.uuid4()
        self.worker_id = uuid.uuid4()
        self.worker_pid = os.getpid()
        
        self.context = zmq.Context()
                
    def tearDown(self):
        #self.test_client.shutdown_workers()
        if self.test_client._monitor_thread is not None and self.test_client._monitor_thread.is_alive():
            self.test_client.shutdown()
            self.test_client._close_signal_sockets()
        del self.context, self.server_id, self.node_id
        
    def check_properly_shut_down(self):
        assert not self.test_client.workers
        assert self.test_client._monitor_thread is not None and not self.test_client._monitor_thread.is_alive()
        assert self.test_client._taskfwd_thread is not None and not self.test_client._taskfwd_thread.is_alive()
        assert self.test_client._rslfwd_thread is not None and not self.test_client._rslfwd_thread.is_alive()
        
    def test_spawn_workers(self):
        '''Client: spawns workers'''
        self.test_client.startup()
        sockdelay()
        
        try:
            for (_pid,proc) in self.test_client.workers.items():
                assert proc.is_alive()
        finally:        
            self.test_client.shutdown()
        
    def test_shutdown_workers(self):
        '''Client: shuts down workers'''
        self.test_client.startup()
        sockdelay()
        self.test_client.shutdown()
        sockdelay()
        assert not self.test_client.workers

    
    def test_startup(self):
        '''Client: all threads start up'''
        self.test_client.startup()
        
        try:
            assert self.test_client.workers
            assert self.test_client._taskfwd_thread is not None and self.test_client._taskfwd_thread.is_alive()
            assert self.test_client._rslfwd_thread is not None and self.test_client._rslfwd_thread.is_alive()
            assert self.test_client._monitor_thread is not None and self.test_client._monitor_thread.is_alive()
        finally:
            self.test_client.shutdown()
                 
    def test_shutdown_internal(self):
        '''Client: internal shutdown signal causes shutdown'''
        self.test_client.startup()
        sockdelay()
        self.test_client.shutdown()
        self.test_client._wait_for_shutdown()
        self.check_properly_shut_down()
            
    def test_shutdown_external(self):
        '''Client: external shutdown signal causes shutdown'''
        
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
        self.check_properly_shut_down()

    def test_missing_server_shutdown(self):
        '''Client: missing server causes shutdown'''
        announce_socket = self.context.socket(zmq.PUB)
        announce_socket.bind(self.ann_endpoint)
        
        self.test_client.server_heartbeat_interval = 0.1
        self.test_client.startup()
        sockdelay()
        announce_socket.send(zwm.MSG_PING) # to start the count
        sockdelay()
        announce_socket.close(linger=0)
        self.test_client._wait_for_shutdown()
        self.check_properly_shut_down()
        
    @timed(2)
    def test_task_forward_lowlevel(self):
        '''Client: forwards task requests upstream'''
        
        self.test_client.startup()
        self.test_client._shutdown_all_workers() # so they won't interfere with this test
        
        worker_task_socket = self.context.socket(zmq.REQ)
        worker_task_socket.connect(self.test_client.worker_task_endpoint)
        
        upstream_task_socket = self.context.socket(zmq.REP)
        upstream_task_socket.bind(self.test_client.upstream_task_endpoint)
        
        task_id = uuid.uuid4()
        
        try:
            worker_task_socket.send_pyobj((zwm.MSG_TASK_REQUEST, self.node_id, self.worker_id, self.worker_pid), flags=zmq.SNDMORE)
            worker_task_socket.send_pyobj((zwm.MSG_TASK_REQUEST, self.server_id, self.worker_id, None))
            message = upstream_task_socket.recv_pyobj()
            
            # did the message arrive appropriately?
            assert message == (zwm.MSG_TASK_REQUEST, self.server_id, self.worker_id, None)
            
            # send response
            upstream_task_socket.send_pyobj((zwm.MSG_TASK_AVAILABLE, self.server_id, self.worker_id, task_id), flags=zmq.SNDMORE)
            upstream_task_socket.send_pyobj((identity, (1,), {}))
            
            # read response
            frames = worker_task_socket.recv_multipart()
            node_header = zwm.unpickle_frame(frames[0])
            task_header = zwm.unpickle_frame(frames[1])
            task_data = zwm.unpickle_frame(frames[2])
            
            assert node_header == (zwm.MSG_TASK_AVAILABLE, self.node_id, self.worker_id, self.worker_pid)
            assert task_header == (zwm.MSG_TASK_AVAILABLE, self.server_id, self.worker_id, task_id)
            assert task_data == (identity, (1,), {})
                
        finally:
            worker_task_socket.close(linger=0)
            upstream_task_socket.close(linger=0)

    @timed(2)
    def test_result_forward_lowlevel(self):
        '''Client: forwards results upstream'''
        
        self.test_client.startup()
        self.test_client._shutdown_all_workers() # so they won't interfere with this test
        
        worker_result_socket = self.context.socket(zmq.REQ)
        worker_result_socket.connect(self.test_client.worker_result_endpoint)
        
        upstream_result_socket = self.context.socket(zmq.REP)
        upstream_result_socket.bind(self.test_client.upstream_result_endpoint)
        
        task_id = uuid.uuid4()
        
        # Put a fake task indicator on our active tasks list
        self.test_client.worker_active_tasks[self.worker_pid] = task_id

        try:
            worker_result_socket.send_pyobj((zwm.MSG_RESULT_SUBMISSION, self.node_id, self.worker_id, self.worker_pid, task_id),
                                            flags=zmq.SNDMORE)
            worker_result_socket.send_pyobj((zwm.MSG_RESULT_SUBMISSION, self.server_id, self.worker_id, task_id),
                                            flags=zmq.SNDMORE)
            worker_result_socket.send_pyobj((zwm.RESULT_TYPE_RETVAL, 1))
            worker_result_socket.recv()
            
            frames = upstream_result_socket.recv_multipart()
            assert zwm.unpickle_frame(frames[0]) == (zwm.MSG_RESULT_SUBMISSION, self.server_id, self.worker_id, task_id)
            assert zwm.unpickle_frame(frames[1]) == (zwm.RESULT_TYPE_RETVAL, 1)
            upstream_result_socket.send('')
            
                
        finally:
            worker_result_socket.close(linger=0)
            upstream_result_socket.close(linger=0)
            
    @timed(2)
    def test_autokill_hung_worker(self):
        '''Client: hung worker is terminated'''
        upstream_task_socket = self.context.socket(zmq.REP)
        upstream_task_socket.bind(self.task_endpoint)
        upstream_result_socket = self.context.socket(zmq.REP)
        upstream_result_socket.bind(self.result_endpoint)
        
        self.test_client.worker_task_timeout = 0.2
        self.test_client.hangcheck = 0.05
        self.test_client.startup()
        
        task_id = uuid.uuid4()

        try:
            # receive task request
            task_req_message = upstream_task_socket.recv_pyobj()
            (_, server_id, worker_id, _) = task_req_message
            
            # dispatch task        
            upstream_task_socket.send_pyobj((zwm.MSG_TASK_AVAILABLE, server_id, worker_id, task_id), flags=zmq.SNDMORE)
            upstream_task_socket.send_pyobj((will_busyhang_uninterruptible, (), {}))
            
            # Receive result
            frames = upstream_result_socket.recv_multipart()
            upstream_result_socket.send('')
            header = zwm.unpickle_frame(frames[0])
            (result_type, payload) = zwm.unpickle_frame(frames[1])
            
            assert header[0] == zwm.MSG_RESULT_SUBMISSION
            assert result_type == zwm.RESULT_TYPE_EXCEPTION
            assert isinstance(payload[0], zwm.WorkerTerminated)
            
        finally:
            upstream_task_socket.close(linger=0)
            upstream_result_socket.close(linger=0)
        
            
class TestCoordinated(CommonParallelTests,CommonWorkManagerTests):
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
        self.test_client._close_signal_sockets()
        self.test_master._close_signal_sockets()
    
    @timed(2)
    def test_sigint_shutdown(self):
        '''Coordination: SIGINT shuts down client and server'''
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
    
    @timed(5)
    def test_stress(self):
        '''Coordination: many small tasks don't crash or lock'''
        N = 1024
        futures = self.work_manager.submit_many([(busy_identity, (n,), {}) for n in xrange(N)])
        self.work_manager.wait_all(futures)
        results = set(future.get_result() for future in futures)
        assert results == set(xrange(N))

           
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

    def test_environ_empty(self):
        with ZMQWorkManager.from_environ() as work_manager:
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
                    
    def test_worker_ids(self):
        work_manager = ZMQWorkManager()
        with work_manager:
            futures = work_manager.submit_many([(get_process_index, (), {})] * work_manager.n_workers)
            work_manager.wait_all(futures)
            results = set(future.get_result() for future in futures)
            assert results == set(str(n) for n in xrange(work_manager.n_workers))
            
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
