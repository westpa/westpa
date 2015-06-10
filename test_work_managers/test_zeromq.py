'''
Created on Jun 2, 2015

@author: mzwier
'''

from __future__ import division, print_function; __metaclass__ = type
import os, signal, tempfile, time, sys, multiprocessing, uuid, socket, json, re
import cPickle as pickle

from work_managers import WMFuture
from work_managers import zeromq as zwm
from work_managers.zeromq import ZMQMaster, ZMQWorker, ZMQWorkManager
from work_managers.zeromq.core import Message, Task, Result, ZMQCore
from tsupport import *

from contextlib import contextmanager

import zmq

import nose.tools
from nose.tools import raises, nottest, timed, assert_raises #@UnresolvedImport
from unittest import skip

import logging
logging.basicConfig(level=logging.DEBUG)

SETUP_WAIT = 0.010

# Amount of time to wait prior to executing tearDown(), to ensure that shutdown
# message don't get lost.
# The original value here (0.010 s = 10 ms) is probably quite generous. 
TEARDOWN_WAIT = 0.010

# How long to wait to let shutdown signals sort themselves out
SHUTDOWN_WAIT = 1

def sockdelay():
    '''Delay for slightly longer than the default auto-reconnect time for ZeroMQ (100 ms)'''
    time.sleep(0.2)
    
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

class ZMQTestBase(object):
    '''Support routines'''
    
    # default endpoint type for tests whose transport is not otherwise specified
    endpoint_type = 'ipc' 
    
    def make_ipc_endpoint(self):
        endpoint = randipc()
        try:
            self._endpoints.append(endpoint)
        except AttributeError:
            self._endpoints = [endpoint]
        return endpoint
    
    def make_tcp_endpoint(self):
        return 'tcp://127.0.0.1:{}'.format(randport())
    
    def make_endpoint(self):
        try:
            endpoint_type = self.endpoint_type
        except AttributeError:
            endpoint_type = 'ipc'
        
        if endpoint_type == 'ipc':
            return self.make_ipc_endpoint()
        elif endpoint_type == 'tcp':
            return self.make_tcp_endpoint()
        else:
            raise ValueError('invalid endpoint type set')
            
    def cleanup_endpoints(self):
        for endpoint in self._endpoints:
            try:
                os.unlink(endpoint[6:])
            except OSError:
                pass
        del self._endpoints
        
    def setUp(self):
        self._endpoints = []
        self.test_core = ZMQCore()
        self.test_core.context = self.test_context = zmq.Context()
        self.test_core.validation_fail_action = 'raise'
        
    def tearDown(self):
        self.cleanup_endpoints()
        self.test_context.destroy(linger=1)
        del self.test_context
        
    @skip
    def test_meta(self):
        '''Test the testing environment; should shake out the most egregious hangs'''
        pass
    
class TestZMQWorkManager(ZMQTestBase):
    def setUp(self):
        super(TestZMQWorkManager,self).setUp()
        
        self.ann_endpoint = self.make_endpoint()
        self.task_endpoint = self.make_endpoint()
        self.result_endpoint = self.make_endpoint()
        
        
        self.test_wm = ZMQWorkManager(self.task_endpoint, self.result_endpoint, self.ann_endpoint)
        self.test_wm.validation_fail_action = 'raise'
        self.test_wm.startup()
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_wm.signal_shutdown()
        self.test_wm.comm_thread.join()
        
        super(TestZMQWorkManager,self).tearDown()

    @contextmanager
    def task_socket(self):
        socket = self.test_context.socket(zmq.PULL)
        socket.connect(self.test_wm.task_endpoint)
        
        yield socket
        
        socket.close(linger=1)
        
    @contextmanager
    def result_socket(self):
        socket = self.test_context.socket(zmq.PUSH)
        socket.connect(self.test_wm.result_endpoint)
        
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
        '''A shim pretending to be the remote master and sending a task to the WM'''
        self.test_core.send_message(socket, Message.RESULT, result)
            
    # Work manager shuts down on inproc signal
    def test_internal_shutdown(self):
        self.test_wm.signal_shutdown()
        self.test_wm.join()
        assert not self.test_wm.comm_thread.is_alive()
    
    # Work manager sends shutdown announcement downstream
    def test_shutdown_sends_announcement(self):
        subsocket = self.test_context.socket(zmq.SUB)
        subsocket.setsockopt(zmq.SUBSCRIBE,'')
        subsocket.connect(self.test_wm.ann_endpoint)

        time.sleep(SETUP_WAIT)        
        self.test_wm.signal_shutdown()
        time.sleep(SETUP_WAIT)
        
        msgs = self.test_core.recv_all(subsocket)
        assert Message.SHUTDOWN in [msg.message for msg in msgs]
        
    
    # Direct, without pulling ZMQMaster into the fray
    def test_submission_sends_task(self):
        with self.task_socket() as s:
            _future = self.test_wm.submit(identity, (1,))
            assert self.recv_task(s).args == (1,)
    
    # Direct, without pulling ZMQMaster into the fray
    def test_result_returns(self):
        with self.task_socket() as t, self.result_socket() as r:
            future = self.test_wm.submit(identity, (1,))
            self.send_result(r, self.recv_task(t).execute())
        
        assert future.result == 1
        

class TestZMQMaster(ZMQTestBase):
    
    def setUp(self):
        super(TestZMQMaster,self).setUp()
        
        self.upstream_ann_endpoint = self.make_endpoint()
        self.upstream_task_endpoint = self.make_endpoint()
        self.upstream_result_endpoint = self.make_endpoint()
        
        self.downstream_ann_endpoint = self.make_endpoint()
        self.downstream_rr_endpoint = self.make_endpoint()
        
        self.test_master = ZMQMaster(self.upstream_task_endpoint, self.upstream_result_endpoint, self.upstream_ann_endpoint,
                                     self.downstream_rr_endpoint, self.downstream_ann_endpoint)
        self.test_master.startup()
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_master.signal_shutdown()
        self.test_master.comm_thread.join()
        
        super(TestZMQMaster,self).tearDown()
        

    def test_internal_shutdown(self):
        #'''master shuts down on inproc signal'''
        self.test_master.signal_shutdown()
        self.test_master.join()
        assert not self.test_master.comm_thread.is_alive()
    
    
            
