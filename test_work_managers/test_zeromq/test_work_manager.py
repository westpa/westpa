'''
Created on Jun 2, 2015

@author: mzwier
'''

from __future__ import division, print_function; __metaclass__ = type

import time
from work_managers.zeromq import ZMQWorkManager
from work_managers.zeromq.core import Message, Task, Result
from test_work_managers.tsupport import *

from contextlib import contextmanager

import zmq

import nose.tools
from nose.tools import raises, nottest, timed, assert_raises #@UnresolvedImport
from unittest import skip


from . import SETUP_WAIT, TEARDOWN_WAIT, SHUTDOWN_WAIT, BEACON_PERIOD, BEACON_WAIT
from . import ZMQTestBase

    
class TestZMQWorkManager(ZMQTestBase):
    
    '''Tests for the core task dispersal/retrieval and shutdown operations
    (the parts of the WM that do not require ZMQMaster/ZMQWorker).'''
    def setUp(self):
        super(TestZMQWorkManager,self).setUp()
        
        
        self.test_wm = ZMQWorkManager()

        # Set operation parameters 
        self.test_wm.validation_fail_action = 'raise'
        self.test_wm.master_beacon_period = BEACON_PERIOD
        self.test_wm.task_beacon_period = BEACON_PERIOD

        self.rr_endpoint = self.make_endpoint()
        self.ann_endpoint = self.make_endpoint()
        self.test_wm.rr_endpoints.append(self.rr_endpoint)
        self.test_wm.ann_endpoints.append(self.ann_endpoint)
        self.test_wm.startup()
        
        self.test_core.master_id = self.test_wm.master_id
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_wm.signal_shutdown()
        self.test_wm.comm_thread.join()
        
        super(TestZMQWorkManager,self).tearDown()

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
        subsocket.setsockopt(zmq.SUBSCRIBE,'')
        subsocket.connect(self.ann_endpoint)
        
        time.sleep(SETUP_WAIT)
        yield
        time.sleep(SETUP_WAIT)
        
        msgs = self.test_core.recv_all(subsocket)
        try:
            assert message in [msg.message for msg in msgs]
        finally:
            subsocket.close(linger=1)
            
    # Work manager shuts down on inproc signal
    def test_internal_shutdown(self):
        self.test_wm.signal_shutdown()
        self.test_wm.join()
        assert not self.test_wm.comm_thread.is_alive()
    
    # Work manager sends shutdown announcement downstream
    def test_shutdown_sends_announcement(self):
        with self.expect_announcement(Message.SHUTDOWN):
            self.test_wm.signal_shutdown()
        
    def test_master_beacon(self):
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
        
    
    # Direct, without pulling ZMQMaster into the fray
#     def test_submission_sends_task(self):
#         with self.task_socket() as s:
#             _future = self.test_wm.submit(identity, (1,))
#             assert self.recv_task(s).args == (1,)
#     
#     # Direct, without pulling ZMQMaster into the fray
#     def test_result_returns(self):
#         with self.task_socket() as t, self.result_socket() as r:
#             future = self.test_wm.submit(identity, (1,))
#             self.send_result(r, self.recv_task(t).execute())
#         
#         assert future.result == 1
        


    
    
            
