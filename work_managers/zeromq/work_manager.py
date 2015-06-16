'''
Created on Jun 10, 2015

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, Task, Result, ZMQWMEnvironmentError, ZMQWorkerMissing
from worker import ZMQWorker
from work_managers import WorkManager, WMFuture
import multiprocessing

from core import PassiveMultiTimer

import zmq

from collections import deque

class ZMQWorkManager(ZMQCore,WorkManager):
    def __init__(self, n_local_workers=1):
        ZMQCore.__init__(self)
        WorkManager.__init__(self)
        
        # Endpoints for announcements and request/replies
        # We can use IPC for node-local workers and TCP for remote
        # at the same time
        self.ann_endpoints = []
        self.rr_endpoints = []
        
        # Node-local workers (one thread/process each)
        if n_local_workers > 0:
            local_ann_endpoint = self.make_internal_endpoint()
            local_rr_endpoint = self.make_internal_endpoint()
            self.ann_endpoints.append(local_ann_endpoint)
            self.rr_endpoints.append(local_rr_endpoint)
            
            self.local_workers = [ZMQWorker(local_rr_endpoint, local_ann_endpoint) for _n in xrange(n_local_workers)]
            
        else:
            self.local_workers = []
            
        self.local_worker_processes = [multiprocessing.Process(target = worker.startup) for worker in self.local_workers]    
                
        # Futures indexed by task ID
        self.futures = dict()
        
        # Tasks pending distribution
        self.outgoing_tasks = deque()
        
        # Tasks being processed by workers (indexed by worker_id)
        self.assigned_tasks = dict()
        
        # Identity information and last contact from workers
        self.worker_information = dict() # indexed by worker_id
        self.worker_timeouts = PassiveMultiTimer() # indexed by worker_id
        
        # How may seconds between task available announcements
        self.task_beacon_period = 1.0
        
        # Number of seconds between checks to see which workers have timed out
        self.worker_timeout_check = 5.0
        
        # Amount of time to wait for stray requests to arrive so that workers shut down properly
        self.shutdown_timeout = 0.5
        
        self.master_id = self.node_id
        
    def startup(self):
        for process in self.local_worker_processes:
            process.start()
        super(ZMQWorkManager,self).startup()
        
    def submit(self, fn, args=None, kwargs=None):
        future = WMFuture()
        task = Task(fn, args or (), kwargs or {}, task_id = future.task_id)
        self.futures[task.task_id] = future
        self.outgoing_tasks.append(task)
        return future

    def submit_many(self, tasks):
        futures = []        
        for (fn,args,kwargs) in tasks:
            future = WMFuture()
            task = Task(fn, args, kwargs, task_id = future.task_id)
            self.futures[task.task_id] = future
            self.outgoing_tasks.append(task)
            futures.append(future)        
        return futures

    def send_message(self, socket, message, payload=None, flags=0):
        message = Message(message, payload)
        message.master_id = self.node_id
        super(ZMQWorkManager,self).send_message(socket, message, payload, flags)
        
    def handle_result(self, socket, msg):
        self.send_ack(socket,msg)
        with self.message_validation(msg):
            assert msg.message == Message.RESULT
            assert isinstance(msg.payload, Result)
            assert msg.payload.task_id in self.futures
            assert self.assigned_tasks[msg.src_id].task_id == msg.payload.task_id
                        
        result = msg.payload
        
        future = self.futures.pop(result.task_id)
        del self.assigned_tasks[msg.src_id]
        if result.exception is not None:
            future._set_exception(result.exception, result.traceback)
        else:
            future._set_result(result.result)
            
    def handle_task_request(self, socket, msg):
        if not self.outgoing_tasks:
            # No tasks available
            self.send_nak(socket,msg)
        else:
            task = self.outgoing_tasks.popleft()
            
            worker_id = msg.src_id
            self.assigned_tasks[worker_id] = task
            
            self.send_message(socket, Message.TASK, task)
            
    def update_worker_information(self, msg):
        try:
            self.worker_timeouts.reset(msg.src_id)
        except KeyError:
            self.worker_timeouts.add_timer(msg.src_id,self.worker_beacon_period*self.timeout_factor)
        
        if msg.message == Message.IDENTIFY:
            with self.message_validation(msg):
                assert isinstance(msg.payload, dict)
            self.worker_information[msg.src_id] = msg.payload
        
    def check_workers(self):
        expired_worker_ids = self.worker_timeouts.which_expired()
        for expired_worker_id in expired_worker_ids:
            try:
                worker_description = '{!s} ({!s})'.format(expired_worker_id, 
                                                          self.worker_information[expired_worker_id]['description'])
            except KeyError:
                worker_description = str(expired_worker_id)
            
            self.log.error('no contact from worker {}'.format(expired_worker_id, worker_description))
               
            self.remove_worker(expired_worker_id)
                                        
    def remove_worker(self, worker_id):
        try:
            expired_task = self.assigned_tasks.pop(worker_id)
        except KeyError:
            pass
        else:
            self.log.error('aborting task {!r} running on expired worker {!s}'
                           .format(expired_task, worker_id))
            future = self.futures.pop(expired_task.task_id)
            future._set_exception(ZMQWorkerMissing('worker running this task disappeared'))

    
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
        timers.add_timer('worker_timeout_check', self.worker_timeout_check)
        timers.add_timer('startup_timeout', self.startup_timeout)
        timers.reset()
        
        peer_found = False
        
        try:
            # Send a master alive message immediately; it will get discarded if necessary
            self.send_message(ann_socket, Message.MASTER_BEACON)
            
            while True:
                # If a timer is already expired, next_expiration_in() will return 0, which
                # zeromq interprets as infinite wait; so instead we select a 1 ms wait in this
                # case.                
                poll_results = dict(poller.poll((timers.next_expiration_in() or 0.001)*1000))
                
                if poll_results and not peer_found:
                    timers.remove_timer('startup_timeout')
                    peer_found = True
                
                if inproc_socket in poll_results:
                    msgs = self.recv_all(inproc_socket,validate=False)
                    # Check for shutdown; do nothing else if shutdown is signalled
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        self.log.debug('shutdown received')
                        break
                    # Check for any other wake-up messages
                    for msg in msgs:
                        if msg.message == Message.TASKS_AVAILABLE:
                            pass
                
                if rr_socket in poll_results:
                    msg = self.recv_message(rr_socket)
                    self.update_worker_information(msg)
                    
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
                if timers.expired('worker_timeout_check'):
                    self.check_workers()
                    timers.reset('worker_timeout_check')
                if not peer_found and timers.expired('startup_timeout'):
                    self.log.error('startup phase elapsed with no contact from workers; shutting down')
                    while self.futures:
                        future = self.futures.popitem()[1]
                        future._set_exception(ZMQWorkerMissing('no workers available'))
                    break
                
            # Post a shutdown message
            self.log.debug('sending shutdown on ann_socket')
            self.send_message(ann_socket, Message.SHUTDOWN)            
            poller.unregister(inproc_socket)
            
            # Clear incoming queue of requests, to let clients exit request/reply states gracefully
            # (clients will still timeout in these states if necessary)
            timers.add_timer('shutdown', self.shutdown_timeout)
            while not timers.expired('shutdown'):
                poll_results = dict(poller.poll(self.shutdown_timeout / 10 * 1000))
                if rr_socket in poll_results:
                    msg = self.recv_message(rr_socket)
                    self.send_nak(rr_socket, msg)

        finally:
            self.context.destroy(linger=1)
            self.context = None
            self.remove_ipc_endpoints()
    
    
        