from __future__ import division, print_function; __metaclass__ = type
import os, signal, tempfile, time, sys, multiprocessing, uuid, socket, json, re
import cPickle as pickle

from work_managers import WMFuture
from work_managers import zeromq as zwm
from work_managers.zeromq import ZMQWorkManager, ZMQClient, recvall, ZMQServer, ZMQRouter
from work_managers.zeromq.client import ZMQWMProcess
from work_managers.zeromq.router import ZMQDevice
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

class TestZMQServer:
    def setUp(self):
        self.test_client_context = zmq.Context()
        ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.test_master = ZMQServer(task_endpoint, result_endpoint, ann_endpoint, listen_endpoint, 10)

    def tearDown(self):
        self.test_master.shutdown()
        self.test_master._close_signal_sockets()
        self.test_client_context.destroy(linger=0)
        del self.test_master, self.test_client_context
    
    def send_task_request(self, socket):
        socket.send_pyobj((zwm.core.MSG_TASK_REQUEST, None, None, None))
        
    def receive_task(self, socket):
        frames = socket.recv_multipart()
        #(tag, server_id, client_id, task_id) = zwm.unpickle_frame(frames[0])
        #(fn, args, kwargs) = zwm.unpickle_frame(frames[1])
        return zwm.unpickle_frame(frames[1])

    def test_startup(self):
        '''Server: all threads start up successfully'''

        work_manager = self.test_master
        work_manager.startup()

        assert work_manager._dispatch_thread and work_manager._dispatch_thread.is_alive(), 'dispatch thread not started up'
        assert work_manager._receive_thread and work_manager._receive_thread.is_alive(), 'receive thread not started up'
        assert work_manager._announce_thread and work_manager._announce_thread.is_alive(), 'announce thread not started up'
        assert work_manager._listen_thread and work_manager._listen_thread.is_alive(), 'listen thread not started up'

        work_manager.shutdown()

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
        assert not work_manager._listen_thread.is_alive(), 'listen thread still alive'
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
            assert ann == zwm.core.MSG_SHUTDOWN
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
            assert tag == zwm.core.MSG_TASK_UNAVAILABLE
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
            assert tag == zwm.core.MSG_TASK_AVAILABLE
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
            result_socket.send_pyobj((zwm.core.MSG_RESULT_SUBMISSION, server_id, client_id, task_id), flags=zmq.SNDMORE)
            result_socket.send_pyobj((zwm.core.RESULT_TYPE_RETVAL, 1))            
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
            result_socket.send_pyobj((zwm.core.MSG_RESULT_SUBMISSION, server_id, client_id, task_id), flags=zmq.SNDMORE)
            result_socket.send_pyobj((zwm.core.RESULT_TYPE_EXCEPTION, (ValueError(42), 'traceback')))
            future.get_result()
            
        finally:        
            task_socket.close(linger=0)
            result_socket.close(linger=0)

    def test_client_update_ping(self):
        '''Server: receives client updates successfully (as pings or on client initialization)'''

        work_manager = self.test_master
        work_manager.startup()
        listen_endpoint = work_manager.master_listen_endpoint
        update_socket1 = self.test_client_context.socket(zmq.PUSH)
        update_socket1.connect(listen_endpoint)
        update_socket2 = self.test_client_context.socket(zmq.PUSH)
        update_socket2.connect(listen_endpoint)

        assert work_manager.n_workers == 0

        try:
            update_socket1.send_pyobj((zwm.core.MSG_PING, None, 'client1', 2), flags=zmq.SNDMORE)
            update_socket1.send_pyobj((None, None), flags=zmq.SNDMORE)
            update_socket1.send_pyobj((None, None))

            sockdelay()

            assert work_manager.n_workers == work_manager.clients['client1'] == 2, 'expected 2 workers, but counted {}'.format(work_manager.n_workers)

            update_socket2.send_pyobj((zwm.core.MSG_PING, None, 'client2', 4), flags=zmq.SNDMORE)
            update_socket2.send_pyobj((None, None), flags=zmq.SNDMORE)
            update_socket2.send_pyobj((None, None))

            sockdelay()

            assert work_manager.clients['client2'] == 4, 'clients dictionary did not update successfully'
            assert work_manager.n_workers == 6, 'expected 6 workers, but counted {}'.format(work_manager.n_workers)

        finally:
            update_socket1.close(linger=0)
            update_socket2.close(linger=0)

    def test_client_udpate_shutdown(self):
        '''Server: receives client shutdown message and updates appropriately'''

        work_manager = self.test_master
        work_manager.startup()
        listen_endpoint = work_manager.master_listen_endpoint
        update_socket = self.test_client_context.socket(zmq.PUSH)
        update_socket.connect(listen_endpoint)

        try:
            work_manager.clients['test_client'] = 4 #add a 'test_client' with 4 workers
            assert work_manager.n_workers == 4, 'expected 4 workers, but counted {}'.format(work_manager.n_workers)

            update_socket.send_pyobj((zwm.core.MSG_SHUTDOWN, None, 'test_client', 4), flags=zmq.SNDMORE)
            update_socket.send_pyobj((None, None), flags=zmq.SNDMORE)
            update_socket.send_pyobj((None, None))

            sockdelay()

            assert work_manager.n_workers == 0, 'expected 0 workers, but counted {}'.format(work_manager.n_workers)

        finally:
            update_socket.close(linger=0)

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
            assert ann == zwm.core.MSG_PING
            work_manager.shutdown()
            
        finally:
            ann_socket.close(linger=0)

class TestZMQDevice:
    def setUp(self):
        self.upstream_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.downstream_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        upstream_startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))

        self.context = zmq.Context()

        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(upstream_startup_ctl_endpoint)

        self.device = ZMQDevice('',self.upstream_endpoint, self.downstream_endpoint, upstream_startup_ctl_endpoint, context=self.context)
        self.device.startup()

        ctlsocket.recv()

        ctlsocket.close()

    def tearDown(self):
        self.device.shutdown()
        
        del self.device

    def test_forwarding(self):
        '''Device: Forwards downstream request upstream, returns upstream response downstream'''

        context = zmq.Context()

        upstream_socket = context.socket(zmq.REP)
        upstream_socket.bind(self.upstream_endpoint)

        downstream_socket = context.socket(zmq.REQ)
        downstream_socket.connect(self.downstream_endpoint)


        try:
            downstream_socket.send('request')
            request = upstream_socket.recv()
            assert request == 'request'

            upstream_socket.send('reply')
            reply = downstream_socket.recv()
            assert reply == 'reply'

        finally:

            upstream_socket.close()
            downstream_socket.close()

class TestZMQRouter: 
    def setUp(self):
        self.upstream_ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.upstream_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.upstream_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.upstream_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        upstream_args = (self.upstream_task_endpoint, self.upstream_result_endpoint,
                         self.upstream_ann_endpoint, self.upstream_listen_endpoint)

        self.downstream_ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.downstream_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.downstream_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.downstream_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        downstream_args = (self.downstream_task_endpoint, self.downstream_result_endpoint,
                           self.downstream_ann_endpoint, self.downstream_listen_endpoint)

        args = upstream_args + downstream_args
        self.test_router = ZMQRouter(*(upstream_args + downstream_args)) 

    def tearDown(self):
        self.test_router.shutdown()

        del self.test_router

    def test_startup(self):
        '''Router: Starts up successfully'''

        router = self.test_router
        router.startup()

        assert router._ann_thread is not None and router._ann_thread.is_alive()
        assert router._task_device is not None and router._task_device.is_alive()
        assert router._result_device is not None and router._result_device.is_alive()

    def check_properly_shutdown(self):

        router = self.test_router

        assert not router._ann_thread.is_alive(), "Router: Moniter thread still alive"
        assert not router._task_device.is_alive(), "Router: Task Device still alive"
        assert not router._result_device.is_alive(), "Router: Result Device still alive"


    def test_shutdown_internal(self):
        '''Router: Internal shutdown signal causes shutdown'''

        router = self.test_router
        router.startup()
        router.shutdown()

        self.check_properly_shutdown()

    def test_shutdown_external(self):
        '''Router: External shutdown causes shutdown'''
        
        context = zmq.Context()

        upstream_socket = context.socket(zmq.PUB)
        upstream_socket.bind(self.upstream_ann_endpoint)

        self.test_router.startup()
        
        upstream_socket.send('shutdown')

        self.test_router._wait_for_shutdown()
        self.check_properly_shutdown()

        upstream_socket.close()
        context.destroy()

    def test_shutdown_announced(self):
        '''Router: Announces shutdown to clients'''
        self.test_router.startup()

        context = zmq.Context()
        downstream_socket = context.socket(zmq.SUB)
        downstream_socket.setsockopt(zmq.SUBSCRIBE, '')
        downstream_socket.connect(self.downstream_ann_endpoint)
        sockdelay()

        self.test_router._shutdown()

        msg = downstream_socket.recv()

        assert msg == 'shutdown'

        downstream_socket.close()
        context.destroy()

    def test_missing_server_shutdown(self):
        '''Router: Missing server causes shutdown'''

        context = zmq.Context()
        upstream_socket = context.socket(zmq.PUB)
        upstream_socket.bind(self.upstream_ann_endpoint)
        sockdelay()

        self.test_router.startup()
        self.test_router.server_heartbeat_interval = 0.1

        upstream_socket.send(zwm.core.MSG_PING) #start clock
        sockdelay()
        upstream_socket.close()
        context.destroy()

        self.test_router._wait_for_shutdown()
        self.check_properly_shutdown()

    def test_task_device(self):
        '''Router: 'Task' device works correctly'''

        self.test_router.startup()

        context = zmq.Context()

        downstream_socket = context.socket(zmq.REQ)
        downstream_socket.connect(self.downstream_task_endpoint)
 
        upstream_socket = context.socket(zmq.REP)
        upstream_socket.bind(self.upstream_task_endpoint)

        try:
            downstream_socket.send_pyobj((zwm.core.MSG_TASK_REQUEST, None, None, None))
            msg, _, _, _ = upstream_socket.recv_pyobj()
            
            assert msg == 'task'

            upstream_socket.send_pyobj((zwm.core.MSG_TASK_AVAILABLE, None, None, None), flags=zmq.SNDMORE)
            upstream_socket.send_pyobj((None, [], {}))

            frames = downstream_socket.recv_multipart(copy=False)
            assert len(frames) == 2
            msg, _, _, _ = zwm.unpickle_frame(frames[0])

            assert msg == 'task_avail'

        finally:
            downstream_socket.close(linger=0)
            upstream_socket.close(linger=0)

            context.destroy()

    def test_result_device(self):
        '''Router: 'Result' device works correctly'''

        self.test_router.startup()

        context = zmq.Context()

        downstream_socket = context.socket(zmq.REQ)
        downstream_socket.connect(self.downstream_result_endpoint)

        upstream_socket = context.socket(zmq.REP)
        upstream_socket.bind(self.upstream_result_endpoint)

        try:
            downstream_socket.send_pyobj((zwm.core.MSG_RESULT_SUBMISSION, None, None, None), flags=zmq.SNDMORE)
            downstream_socket.send_pyobj((zwm.core.RESULT_TYPE_RETVAL, None))

            frames = upstream_socket.recv_multipart(copy=False)

            msg, _, _, _ = zwm.unpickle_frame(frames[0])
            assert msg == 'result'

            retval, _ = zwm.unpickle_frame(frames[1])
            assert retval == 'retval'

            upstream_socket.send_pyobj((zwm.core.MSG_ACK, None, None, None))
            reply, _, _, _ = downstream_socket.recv_pyobj()
            assert reply == 'ok'

        finally:
            downstream_socket.close(linger=0)
            upstream_socket.close(linger=0)

            context.destroy()

    def test_listen_device(self):
        '''Router: 'Listen' device works correctly'''

        self.test_router.startup()

        context = zmq.Context()

        downstream_socket = context.socket(zmq.PUSH)
        downstream_socket.connect(self.downstream_listen_endpoint)

        upstream_socket = context.socket(zmq.PULL)
        upstream_socket.bind(self.upstream_listen_endpoint)

        try:
            downstream_socket.send_pyobj((zwm.core.MSG_PING, None, 'test_client', 2), flags=zmq.SNDMORE)
            downstream_socket.send_pyobj(('a', 'b'))

            frames = upstream_socket.recv_multipart(copy=False)

            msg, _, client_id, n_workers = zwm.unpickle_frame(frames[0])

            assert msg == 'ping'
            assert client_id == 'test_client'
            assert n_workers == 2

            worker_id1, worker_id2 = zwm.unpickle_frame(frames[1])

            assert worker_id1 == 'a'
            assert worker_id2 == 'b'

        finally:
            downstream_socket.close(linger=0)
            upstream_socket.close(linger=0)
            context.destroy()
       
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
        
        socket.send_pyobj((zwm.core.MSG_TASK_AVAILABLE, self.node_id, client_id, client_pid), flags=zmq.SNDMORE)
        socket.send_pyobj((zwm.core.MSG_TASK_AVAILABLE, self.server_id, client_id, task_id), flags=zmq.SNDMORE)
        socket.send_pyobj((fn, args, kwargs))
        
    def receive_result(self, socket=None):
        socket = socket or self.result_socket
        frames = socket.recv_multipart()
        socket.send_pyobj((zwm.core.MSG_ACK, None, None, None))
        
        (node_tag, node_id, client_id, client_pid, task_id) = zwm.unpickle_frame(frames[0])
        (server_tag, server_id, client_id, task_id) = zwm.unpickle_frame(frames[1])
        (result_type, payload) = zwm.unpickle_frame(frames[2])
        return (result_type, payload)
        
    def test_task(self):
        '''Process: requests and processes single successful task'''
        
        self.reply_with_task(identity, (1,), {})
        result_type, payload = self.receive_result()
        
        assert result_type == zwm.core.RESULT_TYPE_RETVAL
        assert payload == 1
        
    def test_exception(self):
        '''Process: requests and processes single failing task'''
        
        self.reply_with_task(will_fail, (), {})
        result_type, payload = self.receive_result()
        
        assert result_type == zwm.core.RESULT_TYPE_EXCEPTION
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
        
class BaseTestZMQClient:    

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
        assert self.test_client._update_thread is not None and not self.test_client._update_thread.is_alive()
        
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
            assert self.test_client._update_thread is not None and self.test_client._update_thread.is_alive()
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
        announce_socket.send(zwm.core.MSG_PING) # to start the count
        sockdelay()
        announce_socket.close(linger=0)
        self.test_client._wait_for_shutdown()
        self.check_properly_shut_down()

    @timed(2)
    def test_update_startup(self):
        '''Client: sends updates upstream when starting up'''

        listener_socket = self.context.socket(zmq.PULL)
        listener_socket.bind(self.test_client.upstream_listen_endpoint)
        sockdelay()

        self.test_client.startup()

        try:
            frames = listener_socket.recv_multipart(copy=False)

            header = zwm.unpickle_frame(frames[0])

            print(header)
            assert header == (zwm.core.MSG_PING, None, self.test_client.instance_id, self.test_client.n_workers)

        finally:
            listener_socket.close(linger=0)

    @timed(2)
    def test_update_shutdown(self):
        '''Client: sends updates upstream when shutting down'''

        listener_socket = self.context.socket(zmq.PULL)
        listener_socket.bind(self.test_client.upstream_listen_endpoint)
        sockdelay()

        self.test_client.startup(spawn_workers=False)

        listener_socket.recv_multipart(copy=False)

        self.test_client.shutdown()

        try:
            message = listener_socket.recv_pyobj()

            print(message)
            assert message == (zwm.core.MSG_SHUTDOWN, self.test_client.instance_id, None, None)

        finally:
            listener_socket.close(linger=0)
        
    @timed(2)
    def test_task_forward_lowlevel(self):
        '''Client: forwards task requests upstream'''
        
        self.test_client.startup(spawn_workers=False)
        
        worker_task_socket = self.context.socket(zmq.REQ)
        worker_task_socket.connect(self.test_client.worker_task_endpoint)
        
        upstream_task_socket = self.context.socket(zmq.REP)
        upstream_task_socket.bind(self.test_client.upstream_task_endpoint)
        
        task_id = uuid.uuid4()
        
        try:
            worker_task_socket.send_pyobj((zwm.core.MSG_TASK_REQUEST, self.node_id, self.worker_id, self.worker_pid), flags=zmq.SNDMORE)
            worker_task_socket.send_pyobj((zwm.core.MSG_TASK_REQUEST, self.server_id, self.worker_id, None))
            message = upstream_task_socket.recv_pyobj()
            
            # did the message arrive appropriately?
            print(message)
            assert message == (zwm.core.MSG_TASK_REQUEST, self.server_id, self.worker_id, None)
            
            
            # send response
            upstream_task_socket.send_pyobj((zwm.core.MSG_TASK_AVAILABLE, self.server_id, self.worker_id, task_id), flags=zmq.SNDMORE)
            upstream_task_socket.send_pyobj((identity, (1,), {}))
            
            # read response
            frames = worker_task_socket.recv_multipart()
            node_header = zwm.unpickle_frame(frames[0])
            task_header = zwm.unpickle_frame(frames[1])
            task_data = zwm.unpickle_frame(frames[2])
            
            assert node_header == (zwm.core.MSG_TASK_AVAILABLE, self.node_id, self.worker_id, self.worker_pid)
            assert task_header == (zwm.core.MSG_TASK_AVAILABLE, self.server_id, self.worker_id, task_id)
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
            worker_result_socket.send_pyobj((zwm.core.MSG_RESULT_SUBMISSION, self.node_id, self.worker_id, self.worker_pid, task_id),
                                            flags=zmq.SNDMORE)
            worker_result_socket.send_pyobj((zwm.core.MSG_RESULT_SUBMISSION, self.server_id, self.worker_id, task_id),
                                            flags=zmq.SNDMORE)
            worker_result_socket.send_pyobj((zwm.core.RESULT_TYPE_RETVAL, 1))
            worker_result_socket.recv()
            
            frames = upstream_result_socket.recv_multipart()
            assert zwm.unpickle_frame(frames[0]) == (zwm.core.MSG_RESULT_SUBMISSION, self.server_id, self.worker_id, task_id)
            assert zwm.unpickle_frame(frames[1]) == (zwm.core.RESULT_TYPE_RETVAL, 1)
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
        self.test_client.hangcheck_interval = 0.05
        self.test_client.startup()
        
        task_id = uuid.uuid4()

        try:
            # receive task request
            task_req_message = upstream_task_socket.recv_pyobj()
            (_, server_id, worker_id, _) = task_req_message
            
            # dispatch task        
            upstream_task_socket.send_pyobj((zwm.core.MSG_TASK_AVAILABLE, server_id, worker_id, task_id), flags=zmq.SNDMORE)
            upstream_task_socket.send_pyobj((will_busyhang_uninterruptible, (), {}))
            
            # Receive result
            frames = upstream_result_socket.recv_multipart()
            upstream_result_socket.send('')
            header = zwm.unpickle_frame(frames[0])
            (result_type, payload) = zwm.unpickle_frame(frames[1])
            
            assert header[0] == zwm.core.MSG_RESULT_SUBMISSION
            assert result_type == zwm.core.RESULT_TYPE_EXCEPTION
            assert isinstance(payload[0], zwm.WorkerTerminated)
            
        finally:
            upstream_task_socket.close(linger=0)
            upstream_result_socket.close(linger=0)

    def test_knows_number_workers(self):
        '''Client: has an accurate record of the number of workers'''

        assert len(self.test_client.workers) == 0

        initial_n_workers = self.test_client.n_workers

        self.test_client.startup()

        sockdelay()

        assert len(self.test_client.workers) == initial_n_workers

        self.test_client._spawn_worker()
        sockdelay()
        assert len(self.test_client.workers) == initial_n_workers + 1

        for worker in self.test_client.workers.values():

            current_workers = len(self.test_client.workers)
            self.test_client._shutdown_worker(worker)
            sockdelay()
            assert len(self.test_client.workers) == current_workers - 1
        
class TestClientTCPComm(BaseTestZMQClient):
    def setUp(self):
        self.ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
                
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint,
                                     self.listen_endpoint, comm_mode='tcp')
        self.node_id = self.test_client.instance_id
        self.server_id = uuid.uuid4()
        self.worker_id = uuid.uuid4()
        self.worker_pid = os.getpid()
        
        self.context = zmq.Context()
                
class TestClientIPCComm(BaseTestZMQClient):
    def setUp(self):
        self.ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
                
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint,
                                     self.listen_endpoint, comm_mode='ipc')
        self.node_id = self.test_client.instance_id
        self.server_id = uuid.uuid4()
        self.worker_id = uuid.uuid4()
        self.worker_pid = os.getpid()
        
        self.context = zmq.Context()

            
class TestCoordinated(CommonParallelTests,CommonWorkManagerTests):
    def setUp(self):
        self.nprocs = 4        
        
        self.ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        
        self.test_master = ZMQWorkManager(self.nprocs, self.task_endpoint, self.result_endpoint,
                                          self.ann_endpoint, self.listen_endpoint)
        self.test_master.startup()
        
        self.test_client = ZMQClient(self.task_endpoint, self.result_endpoint, self.ann_endpoint,
                                     self.listen_endpoint)
        self.test_client.startup()
        
        self.work_manager = self.test_master
        
        self.context = zmq.Context()
    
    def tearDown(self):
        self.test_client.shutdown()
        self.test_master.shutdown()
        self.test_client._close_signal_sockets()
        self.test_master._close_signal_sockets()

    def test_n_workers(self):
        '''Coordination: server has an accurate count of the number of workers'''

        internal_client = self.work_manager.internal_client
        n_workers_expected = self.test_client.n_workers + internal_client.n_workers
        clients = self.work_manager.clients

        assert clients[internal_client.instance_id] and clients[internal_client.instance_id] == internal_client.n_workers
        sockdelay()
        assert clients[self.test_client.instance_id] and clients[self.test_client.instance_id] == self.test_client.n_workers

        assert self.work_manager.n_workers == n_workers_expected, 'expected {} workers but counted {}'.format(self.work_manager.n_workers)
    
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
            assert not work_manager._listen_thread.is_alive(), 'listen thread still alive'
        finally:
            ann_socket.close(linger=0)
    
    @timed(20)
    def test_stress(self):
        '''Coordination: many small tasks don't crash or lock'''
        N = 1024
        futures = self.work_manager.submit_many([(busy_identity, (n,), {}) for n in xrange(N)])
        self.work_manager.wait_all(futures)
        results = set(future.get_result() for future in futures)
        assert results == set(xrange(N))

class TestCoordinatedRouter(CommonParallelTests,CommonWorkManagerTests):
    def setUp(self):
        self.nprocs = 4        
        
        self.up_ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.up_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.up_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.up_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        up_args = [self.up_task_endpoint, self.up_result_endpoint, self.up_ann_endpoint, self.up_listen_endpoint]

        self.down_ann_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.down_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.down_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.down_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        down_args = [self.down_task_endpoint, self.down_result_endpoint, self.down_ann_endpoint, self.down_listen_endpoint]

        all_args = up_args + down_args
        
        self.test_master = ZMQWorkManager(self.nprocs, self.up_task_endpoint, self.up_result_endpoint,
                                          self.up_ann_endpoint, self.up_listen_endpoint)
        self.test_master.startup()
        
        self.test_client = ZMQClient(self.down_task_endpoint, self.down_result_endpoint, self.down_ann_endpoint,
                                     self.down_listen_endpoint)
        self.test_client.startup()

        self.test_router = ZMQRouter(*all_args)
        self.test_router.startup()
        
        self.work_manager = self.test_master
        
        self.context = zmq.Context()
    
    def tearDown(self):
        self.test_client.shutdown()
        self.test_master.shutdown()
        self.test_router.shutdown()
        self.test_client._close_signal_sockets()
        self.test_master._close_signal_sockets()

    def test_n_workers(self):
        '''Coordination (with router): server has an accurate count of the number of workers'''

        internal_client = self.work_manager.internal_client
        n_workers_expected = self.test_client.n_workers + internal_client.n_workers
        clients = self.work_manager.clients

        assert clients[internal_client.instance_id] and clients[internal_client.instance_id] == internal_client.n_workers
        sockdelay()
        assert clients[self.test_client.instance_id] and clients[self.test_client.instance_id] == self.test_client.n_workers

        assert self.work_manager.n_workers == n_workers_expected, 'expected {} workers but counted {}'.format(self.work_manager.n_workers)
 
    
    @timed(20)
    def test_stress(self):
        '''Coordination (with router): many small tasks don't crash or lock'''
        N = 1024
        futures = self.work_manager.submit_many([(busy_identity, (n,), {}) for n in xrange(N)])
        self.work_manager.wait_all(futures)
        results = set(future.get_result() for future in futures)
        assert results == set(xrange(N))



           

