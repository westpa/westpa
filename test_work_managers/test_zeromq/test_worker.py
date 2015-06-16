from __future__ import division, print_function; __metaclass__ = type

import time
from work_managers.zeromq import ZMQWorker
from work_managers.zeromq.core import Message, Task, Result, TIMEOUT_MASTER_BEACON
from test_work_managers.tsupport import *

from contextlib import contextmanager

import zmq

import nose.tools
from nose.tools import raises, nottest, timed, assert_raises #@UnresolvedImport
from unittest import skip


from . import SETUP_WAIT, TEARDOWN_WAIT, SHUTDOWN_WAIT, BEACON_PERIOD, BEACON_WAIT
from . import ZMQTestBase

class TestZMQWorkerBasic(ZMQTestBase):
    
    #endpoint_type = 'tcp'
    
    '''Tests for the core task dispersal/retrieval and shutdown operations
    (the parts of the WM that do not require ZMQWorker).'''
    def setUp(self):
        super(TestZMQWorkerBasic,self).setUp()

        self.rr_endpoint = self.make_endpoint()
        self.ann_endpoint = self.make_endpoint()
        
        # Need to bind ann_socket here in setup, because if we bind it during 
        # tests, messages get lost.
        self.ann_socket = self.test_context.socket(zmq.PUB)
        self.ann_socket.bind(self.ann_endpoint)
        
        # If we're binding ann_socket, we might as well bind rr_socket
        self.rr_socket = self.test_context.socket(zmq.REP)
        self.rr_socket.bind(self.rr_endpoint)
        
        self.test_worker = ZMQWorker(self.rr_endpoint, self.ann_endpoint)
        self.test_worker.validation_fail_action = 'raise'
        self.test_worker.shutdown_timeout = 0.5
        self.test_worker.master_beacon_period = BEACON_PERIOD
        self.test_worker.startup()
        
        self.test_core.master_id = self.test_core.node_id
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_worker.signal_shutdown()
        self.test_worker.comm_thread.join()
        
        super(TestZMQWorkerBasic,self).tearDown()

    def send_task(self, task):
        self.test_core.send_message(self.ann_socket, Message.TASKS_AVAILABLE)
        msg = self.test_core.recv_message(self.rr_socket)
        assert msg.message == Message.TASK_REQUEST
        self.test_core.send_message(self.rr_socket, Message.TASK, payload=task)
        
    def recv_result(self):
        msg = self.test_core.recv_message(self.rr_socket)
        self.test_core.send_ack(self.rr_socket,msg)
        assert msg.message == Message.RESULT
        assert isinstance(msg.payload, Result)
        return msg.payload
         
    def roundtrip_task(self, task):
        self.send_task(task)
        return self.recv_result()


    def test_meta(self):
        pass
    
    def test_executor_alive(self):
        assert self.test_worker.executor_process.is_alive()
        
    def test_executor_shuts_down_immediately(self):
        self.test_worker.shutdown_executor()
        assert not self.test_worker.executor_process.is_alive()
        
    def test_shutdown_on_announcement(self):
        self.test_core.send_message(self.ann_socket, Message.SHUTDOWN)
        self.test_worker.join()
        assert not self.test_worker.executor_process.is_alive()
    
    def test_responds_to_task_avail(self):
        self.test_core.send_message(self.ann_socket, Message.TASKS_AVAILABLE)
        msg = self.test_core.recv_message(self.rr_socket)
        self.test_core.send_nak(self.rr_socket,msg)
        assert msg.message == Message.TASK_REQUEST
        
    def test_shutdown_on_master_disappearance(self):
        self.test_core.send_message(self.ann_socket, Message.RECONFIGURE_TIMEOUT, (TIMEOUT_MASTER_BEACON, 0.01))
        time.sleep(0.02)
        self.test_worker.join()
        assert not self.test_worker.executor_process.is_alive()        
        
    def test_worker_processes_task(self):
        r = random_int()
        task = Task(identity, (r,), {})
        rsl = self.roundtrip_task(task) 
        assert rsl.result == r
                
    def test_worker_processes_exception(self):
        task = Task(will_fail, (), {})
        rsl = self.roundtrip_task(task)
        assert isinstance(rsl.exception, ExceptionForTest)
        
    def test_hung_worker_interruptible(self):
        task = Task(will_busyhang, (), {})
        self.send_task(task)
        time.sleep(1.0)
        self.test_core.send_message(self.ann_socket, Message.SHUTDOWN)
        self.test_worker.join()
    
    def test_hung_worker_uninterruptible(self):
        task = Task(will_busyhang_uninterruptible, (), {})
        self.send_task(task)
        time.sleep(1.0)
        self.test_core.send_message(self.ann_socket, Message.SHUTDOWN)
        self.test_worker.join()
