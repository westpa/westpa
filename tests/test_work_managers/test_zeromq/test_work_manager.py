'''
Created on Jun 2, 2015

@author: mzwier
'''


import time
from work_managers.zeromq import ZMQWorkManager, ZMQWorker, ZMQWorkerMissing
from work_managers.zeromq.core import Message, Task
from test_work_managers.tsupport import *

from contextlib import contextmanager

import zmq

import nose.tools
from nose.tools import raises, nottest, timed, assert_raises #@UnresolvedImport
from unittest import skip


from . import SETUP_WAIT, TEARDOWN_WAIT, SHUTDOWN_WAIT, BEACON_PERIOD, BEACON_WAIT
from . import ZMQTestBase

    
class TestZMQWorkManagerBasic(ZMQTestBase):
    
    '''Tests for the core task dispersal/retrieval and shutdown operations
    (the parts of the WM that do not require ZMQWorker).'''
    def setUp(self):
        super(TestZMQWorkManagerBasic,self).setUp()
        
        
        self.test_wm = ZMQWorkManager(n_local_workers=0)

        # Set operation parameters 
        self.test_wm.validation_fail_action = 'raise'
        self.test_wm.master_beacon_period = BEACON_PERIOD
        self.test_wm.task_beacon_period = BEACON_PERIOD

        self.rr_endpoint = self.make_endpoint()
        self.ann_endpoint = self.make_endpoint()
        self.test_wm.downstream_rr_endpoint = self.rr_endpoint
        self.test_wm.downstream_ann_endpoint = self.ann_endpoint
        self.test_wm.startup()
        
        self.test_core.master_id = self.test_wm.master_id
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_wm.signal_shutdown()
        self.test_wm.comm_thread.join()
        
        super(TestZMQWorkManagerBasic,self).tearDown()

    @contextmanager
    def task_socket(self):
        socket = self.test_context.socket(zmq.PULL)
        socket.connect(self.test_wm.wm_task_endpoint)
        
        yield socket
        
        socket.close(linger=1)
        
    @contextmanager
    def result_socket(self):
        socket = self.test_context.socket(zmq.PUSH)
        socket.connect(self.test_wm.wm_result_endpoint)
        
        yield socket
        
        socket.close(linger=1)

    def recv_task(self, socket):
        '''A shim pretending to receive a task from a remote master'''
        msg = self.test_core.recv_message(socket)
        assert msg.message == Message.TASK
        assert msg.payload is not None
        assert isinstance(msg.payload, Task)
        return msg.payload
    
    def send_result(self, socket, result):
        '''A shim pretending to be the remote worker and sending a task to the WM'''
        self.test_core.send_message(socket, Message.RESULT, result)
        
    @contextmanager
    def rr_socket(self):
        socket = self.test_context.socket(zmq.REQ)
        socket.connect(self.rr_endpoint)
        
        yield socket
        
        socket.close(linger=1)
        
    @contextmanager
    def expect_announcement(self, message):
        subsocket = self.test_context.socket(zmq.SUB)
        subsocket.setsockopt(zmq.SUBSCRIBE,b'')
        subsocket.connect(self.ann_endpoint)
        
        time.sleep(SETUP_WAIT)
        yield
        time.sleep(SETUP_WAIT)
        
        msgs = self.test_core.recv_all(subsocket)
        try:
            assert message in [msg.message for msg in msgs]
        finally:
            subsocket.close(linger=0)
            
    def discard_announcements(self):
        subsocket = self.test_context.socket(zmq.SUB)
        subsocket.setsockopt(zmq.SUBSCRIBE,b'')
        subsocket.connect(self.ann_endpoint)
        
        time.sleep(SETUP_WAIT)
        self.test_core.recv_all(subsocket)
        subsocket.close(linger=0)        
            
    # Work manager shuts down on inproc signal
    def test_internal_shutdown(self):
        self.test_wm.signal_shutdown()
        self.test_wm.join()
        assert not self.test_wm.comm_thread.is_alive()
    
    # Work manager sends shutdown announcement downstream
    def test_shutdown_sends_announcement(self):
        with self.expect_announcement(Message.SHUTDOWN):
            self.test_wm.signal_shutdown()
        
# This won't work, because initial beacon is discarded if no clients are connected
#     def test_immediate_master_beacon(self):
#         with self.expect_announcement(Message.MASTER_BEACON):
#             time.sleep(BEACON_WAIT)
    
    @skip       
    def test_delayed_master_beacon(self):
        self.discard_announcements()
        with self.expect_announcement(Message.MASTER_BEACON):
            time.sleep(BEACON_WAIT)
    
            
    def test_task_avail_beacon(self):
        with self.expect_announcement(Message.TASKS_AVAILABLE):
            self.test_wm.submit(identity, (random_int(),))
            time.sleep(BEACON_WAIT)
            
    def test_task_nak(self):
        with self.rr_socket() as s:
            self.test_core.send_message(s,Message.TASK_REQUEST)
            msg = self.test_core.recv_message(s)
            assert msg.message == Message.NAK
            
    def test_task_send(self):
        r = random_int()
        self.test_wm.submit(identity, (r,))
        with self.rr_socket() as s:
            self.test_core.send_message(s,Message.TASK_REQUEST)
            msg = self.test_core.recv_message(s)
            assert msg.message == Message.TASK
            assert isinstance(msg.payload, Task)
            assert msg.payload.args == (r,)
            
    def test_result_return(self):
        r = random_int()
        future = self.test_wm.submit(identity, (r,))
        with self.rr_socket() as s:
            self.test_core.send_message(s,Message.TASK_REQUEST)
            msg = self.test_core.recv_message(s)
            task = msg.payload
            result = task.execute()
            self.test_core.send_message(s, Message.RESULT, result)
        assert future.result == r

class BaseInternal(ZMQTestBase,CommonWorkManagerTests):
    def setUp(self):
        super(BaseInternal,self).setUp()
        
        
        self.test_wm = ZMQWorkManager(n_local_workers=self.n_workers)
        for worker in self.test_wm.local_workers:
            worker.validation_fail_action = 'raise'
            worker.shutdown_timeout = 0.5

        # Set operation parameters 
        self.test_wm.validation_fail_action = 'raise'
        self.test_wm.master_beacon_period = BEACON_PERIOD
        self.test_wm.task_beacon_period = BEACON_PERIOD
        self.test_wm.startup_timeout = 1.0
        self.test_wm.shutdown_timeout = 0.5

        self.test_wm.startup()
        
        self.test_core.master_id = self.test_wm.master_id
        
        self.work_manager = self.test_wm
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_wm.signal_shutdown()
        self.test_wm.comm_thread.join()
        
        super(BaseInternal,self).tearDown()

    
    def test_worker_startup(self):
        time.sleep(0.1)
        for worker_process in self.test_wm.local_worker_processes:
            assert worker_process.is_alive()
        
    def test_processes_task(self):
        r = random_int()
        future = self.test_wm.submit(identity,(r,),{})
        assert future.get_result() == r
        
    def test_two_tasks(self):
        r1 = random_int()
        r2 = random_int()
        f1 = self.test_wm.submit(identity,(r1,),{})
        f2 = self.test_wm.submit(identity,(r2,),{})
        assert f1.result == r1
        assert f2.result == r2
    
    def test_processes_exception(self):
        future = self.test_wm.submit(will_fail, (), {})
        assert isinstance(future.get_exception(), ExceptionForTest)
        
    def test_hung_worker_interruptible(self):
        self.test_wm.submit(will_busyhang, (), {})
        time.sleep(1.0)
        
    def test_hung_worker_uninterruptible(self):
        self.test_wm.submit(will_busyhang_uninterruptible, (), {})
        time.sleep(1.0)
        
class TestZMQWorkManagerInternalNone(ZMQTestBase):
    n_workers = 0
 
    def setUp(self):
        super(TestZMQWorkManagerInternalNone,self).setUp()
        
        
        self.test_wm = ZMQWorkManager(n_local_workers=self.n_workers)
        for worker in self.test_wm.local_workers:
            worker.validation_fail_action = 'raise'
            worker.shutdown_timeout = 0.5

        # Set operation parameters 
        self.test_wm.validation_fail_action = 'raise'
        self.test_wm.master_beacon_period = BEACON_PERIOD
        self.test_wm.task_beacon_period = BEACON_PERIOD
        self.test_wm.startup_timeout = 1.0
        self.test_wm.shutdown_timeout = 0.5

        self.test_wm.startup()
        
        self.test_core.master_id = self.test_wm.master_id
        
        self.work_manager = self.test_wm
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_wm.signal_shutdown()
        self.test_wm.comm_thread.join()
        
        super(TestZMQWorkManagerInternalNone,self).tearDown()
 
    
    def test_shutdown_without_workers(self):
        time.sleep(1.5)
        assert not self.test_wm.comm_thread.is_alive()
        
    def test_shutdown_without_workers_after_submission(self):
        self.test_wm.submit(identity, (1,), {})
        time.sleep(1.5)
        assert not self.test_wm.comm_thread.is_alive()
        
    def test_shutdown_without_workers_raises_future_error(self):
        future = self.test_wm.submit(identity, (1,), {})
        time.sleep(1.5)
        assert isinstance(future.get_exception(),ZMQWorkerMissing)
                 
class TestZMQWorkManagerInternalSingle(BaseInternal):
    n_workers = 1
    
class TestZMQWorkManagerInternalMultiple(BaseInternal):
    n_workers = 4

class BaseExternal(ZMQTestBase,CommonWorkManagerTests):

    def setUp(self):
        super(BaseExternal,self).setUp()
        
        
        self.test_wm = ZMQWorkManager(n_local_workers=0)
        
        # Set operation parameters 
        self.test_wm.validation_fail_action = 'raise'
        self.test_wm.master_beacon_period = BEACON_PERIOD
        self.test_wm.task_beacon_period = BEACON_PERIOD
        self.test_wm.worker_beacon_period = BEACON_PERIOD
        self.test_wm.worker_timeout_check = BEACON_WAIT

        self.rr_endpoint = self.make_endpoint()
        self.ann_endpoint = self.make_endpoint()
        self.test_wm.downstream_rr_endpoint = self.rr_endpoint
        self.test_wm.downstream_ann_endpoint = self.ann_endpoint
        self.test_wm.startup()

        self.workers = [ZMQWorker(self.rr_endpoint, self.ann_endpoint)]        
        for worker in self.workers:
            worker.validation_fail_action = 'raise'
            worker.master_beacon_period = BEACON_PERIOD
            worker.shutdown_timeout = 0.5
            worker.startup()
        
        
        self.test_core.master_id = self.test_wm.master_id
        
        self.work_manager = self.test_wm
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_wm.signal_shutdown()
        self.test_wm.comm_thread.join()
        
        super(BaseExternal,self).tearDown()
    
class TestZMQWorkManagerExternalSingle(BaseExternal):
    n_workers = 1
    
class TestZMQWorkManagerExternalMultiple(BaseExternal):
    n_workers = 4

    
