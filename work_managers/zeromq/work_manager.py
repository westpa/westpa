'''
Created on Jun 10, 2015

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, Task, Result, ZMQWMEnvironmentError, ZMQWorkerMissing
from work_managers import WorkManager, WMFuture

import zmq

class ZMQWorkManager(ZMQCore,WorkManager):
    def __init__(self, task_endpoint, result_endpoint, ann_endpoint):
        ZMQCore.__init__(self)
        WorkManager.__init__(self)
        
        # ID of ZMQMaster (usually running in separate thread/process)
        self.master_id = None
        self.task_endpoint = task_endpoint
        self.result_endpoint = result_endpoint
        self.ann_endpoint = ann_endpoint
        
        # Futures indexed by task ID
        self.futures = dict()
        
    def submit(self, fn, args=None, kwargs=None):
        future = WMFuture()
        task = Task(fn, args or (), kwargs or {}, task_id = future.task_id)
        self.futures[task.task_id] = future
        self.send_message(self.task_socket, Message.TASK, payload=task)
        return future
                

    def send_shutdown_message(self, signal=None):
        self.send_message(self.ann_socket, Message.SHUTDOWN)
        
    
    def handle_result_message(self, msg):
        with self.message_validation(msg):
            assert msg.message == Message.RESULT
            assert isinstance(msg.payload, Result)
            assert msg.payload.task_id in self.futures
                        
        result = msg.payload
        
        future = self.futures.pop(result.task_id)
        if result.exception is not None:
            future._set_exception(result.exception, result.traceback)
        else:
            future._set_result(result.result)
    
    def comm_loop(self):
        self.context = zmq.Context()
        
        self.task_socket = self.context.socket(zmq.PUSH)
        self.result_socket = self.context.socket(zmq.PULL)
        self.ann_socket = self.context.socket(zmq.PUB)
        
        self.task_socket.bind(self.task_endpoint)
        self.result_socket.bind(self.result_endpoint)
        self.ann_socket.bind(self.ann_endpoint)

        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
        inproc_socket.bind(self.inproc_endpoint)
        
        poller = zmq.Poller()
        poller.register(inproc_socket, zmq.POLLIN)
        poller.register(self.result_socket, zmq.POLLIN)
        
        try:
            while True:
                poll_results = dict(poller.poll())
                
                # Check for shutdown
                if inproc_socket in poll_results:
                    msgs = self.recv_all(inproc_socket,validate=False)
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        return
                
                # Process results
                if self.result_socket in poll_results:
                    for msg in self.recv_all(self.result_socket):
                        self.handle_result_message(msg)
                    
            
        finally:
            self.send_shutdown_message()
            self.context.destroy(linger=1)
    
    
        