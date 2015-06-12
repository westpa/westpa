'''
Created on Jun 10, 2015

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, Task, Result, ZMQWMEnvironmentError, ZMQWorkerMissing
from master import ZMQMaster
from work_managers import WorkManager, WMFuture

from core import PassiveMultiTimer

import zmq

from collections import deque

class ZMQWorkManager(ZMQCore,WorkManager):
    def __init__(self, n_workers=1):
        ZMQCore.__init__(self)
        WorkManager.__init__(self)
        
        # Endpoints for announcements and request/replies
        # We can use IPC for node-local workers and TCP for remote
        # at the same time
        self.ann_endpoints = []
        self.rr_endpoints = []
        
        # Node-local workers (one thread/process each)
        self.local_workers = []
                
        # Futures indexed by task ID
        self.futures = dict()
        
        # Tasks pending distribution
        self.outgoing_tasks = deque()
        
        # Tasks being processed by workers (indexed by worker_id)
        self.assigned_tasks = dict()        
        
        # How may seconds between task available announcements
        self.task_beacon_period = 1.0
        
        # Number of seconds between checks to see which workers have timed out
        self.worker_timeout_check = 5.0
        
        self.master_id = self.node_id
        
    def submit(self, fn, args=None, kwargs=None):
        future = WMFuture()
        task = Task(fn, args or (), kwargs or {}, task_id = future.task_id)
        self.futures[task.task_id] = future
        self.outgoing_tasks.append(task)
        self.send_inproc_message(Message.TASKS_AVAILABLE)
        return future

    def submit_many(self, tasks):
        futures = []
        task_socket = self.context.socket(zmq.PUSH)
        task_socket.connect(self.wm_task_endpoint)
        
        for (fn,args,kwargs) in tasks:
            future = WMFuture()
            task = Task(fn, args, kwargs, task_id = future.task_id)
            self.futures[task.task_id] = future
            self.outgoing_tasks.append(task)
            futures.append(future) 
        self.send_inproc_message(Message.TASKS_AVAILABLE)       
        return futures

    def send_message(self, socket, message, payload=None, flags=0):
        message = Message(message, payload)
        message.master_id = self.node_id
        super(ZMQWorkManager,self).send_message(socket, message, payload, flags)
        
    def handle_result(self, socket, msg):
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
            
        self.send_ack(socket,msg)
        
    def handle_task_request(self, socket, msg):
        if not self.outgoing_tasks:
            # No tasks available
            self.send_nak(socket,msg)
        else:
            task = self.outgoing_tasks.popleft()
            
            worker_id = msg.src_id
            self.assigned_tasks[worker_id] = task
            
            self.send_message(socket, Message.TASK, task)
    
    def comm_loop(self):
        self.context = zmq.Context()
        
        rr_socket = self.context.socket(zmq.REP)
        ann_socket = self.context.socket(zmq.PUB)
        
        for endpoint in self.rr_endpoints:
            rr_socket.bind(endpoint)
            
        for endpoint in self.ann_endpoints:
            ann_socket.bind(endpoint)

        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
        inproc_socket.bind(self.inproc_endpoint)
        
        poller = zmq.Poller()
        poller.register(inproc_socket, zmq.POLLIN)
        poller.register(rr_socket, zmq.POLLIN)
        
        timers = PassiveMultiTimer()
        timers.add_timer('tasks_avail', self.task_beacon_period)
        timers.add_timer('master_beacon', self.master_beacon_period)
        timers.reset()
        
        try:
            # Send a master alive message immediately
            self.send_message(ann_socket, Message.MASTER_BEACON)
            
            while True:
                # If a timer is already expired, next_expiration_in() will return 0, which
                # zeromq interprets as infinite wait; so instead we select a 1 ms wait in this
                # case.                
                poll_results = dict(poller.poll((timers.next_expiration_in() or 0.001)*1000))
                
                if inproc_socket in poll_results:
                    msgs = self.recv_all(inproc_socket,validate=False)
                    # Check for shutdown; do nothing else if shutdown is signalled
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        return
                    # Check for any other wake-up messages
                    for msg in msgs:
                        if msg.message == Message.TASKS_AVAILABLE:
                            pass
                
                if rr_socket in poll_results:
                    msg = self.recv_message(rr_socket)
                    
                    if msg.message == Message.TASK_REQUEST:
                        self.handle_task_request(rr_socket, msg)
                    elif msg.message == Message.RESULT:
                        self.handle_result(rr_socket, msg)
                    else:
                        self.send_ack(rr_socket, msg)
                
                if timers.expired('tasks_avail'):
                    if self.outgoing_tasks:
                        self.send_message(ann_socket, Message.TASKS_AVAILABLE)
                    timers.reset('tasks_avail')
                if timers.expired('master_beacon'):
                    self.send_message(ann_socket, Message.MASTER_BEACON)
                    timers.reset('master_beacon')
                
                
        finally:
            self.send_message(ann_socket, Message.SHUTDOWN)
            self.context.destroy(linger=1)
    
    
        