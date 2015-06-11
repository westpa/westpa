'''
Created on Jun 10, 2015

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, Task, Result, ZMQWMEnvironmentError, ZMQWorkerMissing
from master import ZMQMaster
from work_managers import WorkManager, WMFuture


import zmq

class ZMQWorkManager(ZMQCore,WorkManager):
    def __init__(self, n_workers=1):
        ZMQCore.__init__(self)
        WorkManager.__init__(self)

        
        # ID of ZMQMaster (usually running in separate thread/process)
        self.master_id = None
        
        # Endpoints that the work manager uses to communicate with the master
        self.wm_task_endpoint = self.make_internal_endpoint()
        self.wm_result_endpoint = self.make_internal_endpoint()
        self.wm_ann_endpoint = self.make_internal_endpoint()
        
        # ZMQMaster instance
        # Currently, this is always internal to the work manager
        # In theory, it can be run in a separate process
        self.master = ZMQMaster(self.wm_task_endpoint, self.wm_result_endpoint, self.wm_ann_endpoint)                 
        
        
        # Futures indexed by task ID
        self.futures = dict()
        
    def submit(self, fn, args=None, kwargs=None):
        future = WMFuture()
        task = Task(fn, args or (), kwargs or {}, task_id = future.task_id)
        self.futures[task.task_id] = future
        task_socket = self.context.socket(zmq.PUSH)
        task_socket.connect(self.wm_task_endpoint)
        self.send_message(task_socket, Message.TASK, payload=task)
        task_socket.close(linger=1)
        return future

    def submit_many(self, tasks):
        '''Submit a set of tasks to the work manager, returning a list of `WMFuture` objects representing
        pending results. Each entry in ``tasks`` should be a triple (fn, args, kwargs), which will result in
        fn(*args, **kwargs) being executed by a worker. The function ``fn`` and all arguments must be
        picklable; note particularly that off-path modules are not picklable unless pre-loaded in the worker
        process.'''

        futures = []
        task_socket = self.context.socket(zmq.PUSH)
        task_socket.connect(self.wm_task_endpoint)
        
        for (fn,args,kwargs) in tasks:
            future = WMFuture()
            task = Task(fn, args, kwargs, task_id = future.task_id)
            self.futures[task.task_id] = future
            futures.append(future)
            self.send_message(task_socket, Message.TASK, payload=task)
        task_socket.close(linger=1)            
        return futures

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
        
        result_socket = self.context.socket(zmq.PULL)
        ann_socket = self.context.socket(zmq.PUB)
        
        result_socket.bind(self.wm_result_endpoint)
        ann_socket.bind(self.wm_ann_endpoint)

        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
        inproc_socket.bind(self.inproc_endpoint)
        
        poller = zmq.Poller()
        poller.register(inproc_socket, zmq.POLLIN)
        poller.register(result_socket, zmq.POLLIN)
        
        try:
            while True:
                poll_results = dict(poller.poll())
                
                # Check for shutdown
                if inproc_socket in poll_results:
                    msgs = self.recv_all(inproc_socket,validate=False)
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        return
                
                # Process results
                if result_socket in poll_results:
                    for msg in self.recv_all(result_socket):
                        self.handle_result_message(msg)
                    
            
        finally:
            self.send_message(ann_socket, Message.SHUTDOWN)
            self.context.destroy(linger=1)
    
    
        