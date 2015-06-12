'''
Created on May 29, 2015

@author: mzwier
'''

import sys

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, ZMQWMEnvironmentError, PassiveMultiTimer, Task, Result, TIMEOUT_MASTER_BEACON
import threading, multiprocessing, os, signal

import zmq

class ZMQWorker(ZMQCore):
    '''This is the outward facing worker component of the ZMQ work manager. This
    forms the interface to the master. This process cannot hang or crash due to an 
    error in tasks it executes, so tasks are isolated in ZMQExecutor, which 
    communicates with ZMQWorker via (what else?) ZeroMQ.'''
    
    def __init__(self, rr_endpoint, ann_endpoint):
        super(ZMQWorker,self).__init__()
        
        # Upstream endpoints
        self.rr_endpoint = rr_endpoint
        self.ann_endpoint = ann_endpoint
                
        # Downstream endpoints
        self.task_endpoint = self.make_internal_endpoint()
        self.result_endpoint = self.make_internal_endpoint()
        
        self.master_id = None
        self.identified = False
        
        # The task currently being processed
        self.pending_task = None
        
        # Executor process
        
        self.shutdown_timeout = 5.0 # Five second wait between shutdown message and SIGINT and SIGINT and SIGKILL
        self.executor_process = None
            
    def update_master_info(self, msg):
        if self.master_id is None:
            self.master_id = msg.master_id
        self.timers.reset(TIMEOUT_MASTER_BEACON)
        
    def identify(self, rr_socket):
        if self.master_id is None or self.identified or self.timers.expired(TIMEOUT_MASTER_BEACON): return
        self.send_message(rr_socket, Message.IDENTIFY, payload=self.get_identification())
        self.recv_ack(rr_socket)
        self.identified = True
        
    def request_task(self, rr_socket, task_socket):
        if self.master_id is None: return
        elif self.pending_task is not None: return
        elif self.timers.expired(TIMEOUT_MASTER_BEACON): return
        else:
            self.send_message(rr_socket, Message.TASK_REQUEST)
            reply = self.recv_message(rr_socket)
            if reply.message == Message.NAK:
                # No task available
                return 
            else:
                with self.message_validation(reply):
                    assert isinstance(reply.payload, Task)
                    task = reply.payload
                self.pending_task = task
                self.send_message(task_socket, Message.TASK, task)                       
            
    def handle_reconfigure_timeout(self, msg, timers):
        with self.message_validation(msg):
            assert msg.payload is not None
            timer, new_period = msg.payload
        timers.change_duration(timer, new_period)
        timers.reset(timer)
        
    def handle_result(self, result_socket, rr_socket):
        msg = self.recv_message(result_socket)
        with self.message_validation(msg):
            assert msg.message == Message.RESULT
            assert isinstance(msg.payload, Result)
            assert msg.payload.task_id == self.pending_task.task_id
        
        msg.src_id = self.node_id
        self.pending_task = None
        self.send_message(rr_socket, msg)
        self.recv_ack(rr_socket)
    
    def comm_loop(self):
        '''Master communication loop for the worker process.'''
        
        rr_socket = self.context.socket(zmq.REQ)

        ann_socket = self.context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')
        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
        
        task_socket = self.context.socket(zmq.PUSH)
        result_socket = self.context.socket(zmq.PULL)
        
        self.log.info('This is {}'.format(self.node_description))
        
        timers = self.timers = PassiveMultiTimer()
        timers.add_timer(TIMEOUT_MASTER_BEACON, self.master_beacon_period)
        timers.add_timer('worker_beacon', self.worker_beacon_period)
        
        try:
            rr_socket.connect(self.rr_endpoint)            
            ann_socket.connect(self.ann_endpoint)
            inproc_socket.bind(self.inproc_endpoint)
            task_socket.connect(self.task_endpoint)
            result_socket.bind(self.result_endpoint)
            
            poller = zmq.Poller()
            poller.register(ann_socket, zmq.POLLIN)
            poller.register(inproc_socket, zmq.POLLIN)
            poller.register(result_socket, zmq.POLLIN)
            
            timers.reset()
            while True:
                # If a timer is already expired, next_expiration_in() will return 0, which
                # zeromq interprets as infinite wait; so instead we select a 1 ms wait in this
                # case.
                poll_results = dict(poller.poll((timers.next_expiration_in() or 0.001)*1000))
                
                announcements = []   
                
                # Check for internal messages first
                if inproc_socket in poll_results:
                    announcements.extend(self.recv_all(inproc_socket))
          
                # Process announcements
                if ann_socket in poll_results:
                    announcements.extend(self.recv_all(ann_socket))
                    
                #announcements = Message.coalesce_announcements(announcements)
                #self.log.debug('received {:d} announcements'.format(len(announcements)))
                    
                # Check for shutdown messages
                if Message.SHUTDOWN in (msg.message for msg in announcements):
                    self.log.debug('received shutdown message')
                    return
                
                # Handle results, so that we clear ourselves of completed tasks
                # before asking for more
                if result_socket in poll_results:
                    self.handle_result(result_socket, rr_socket)
                    # immediately request another task if available
                    if not timers.expired(TIMEOUT_MASTER_BEACON):
                        self.request_task(rr_socket, task_socket)
                
                # Handle any remaining messages                    
                for msg in announcements:
                    if msg.message == Message.MASTER_BEACON:
                        self.update_master_info(msg)
                        self.identify(rr_socket)
                    elif msg.message == Message.RECONFIGURE_TIMEOUT:
                        self.handle_reconfigure_timeout(msg, timers)
                    elif msg.message == Message.TASKS_AVAILABLE:
                        self.update_master_info(msg)
                        self.request_task(rr_socket,task_socket)
                del announcements
                
                            
                # Process timeouts
                #if timers.expired('worker_beacon'):
                    #self.log.debug('worker_beacon timeout')
                    #self.identify(rr_socket)

                if timers.expired(TIMEOUT_MASTER_BEACON):
                    self.log.error('no contact from master; shutting down')
                    return
                            
        finally:
            self.shutdown_executor()
            self.executor_process.join()
            self.context.destroy(linger=1)
            self.context = None
            
    def shutdown_executor(self):
        if self.context is not None:
            try:
                self.log.debug('sending shutdown task to executor')
                task_socket = self.context.socket(zmq.PUSH)
                task_socket.connect(self.task_endpoint)
                self.send_message(task_socket, Message.SHUTDOWN)
                task_socket.close(linger=1)
            except:
                pass
        
        self.executor_process.join(self.shutdown_timeout)
        if self.executor_process.is_alive():            
            self.log.debug('sending SIGINT to worker process {:d}'.format(self.executor_process.pid))
            os.kill(self.executor_process.pid, signal.SIGINT)
            self.executor_process.join(self.shutdown_timeout)
            if self.executor_process.is_alive():
                self.log.warning('sending SIGKILL to worker process {:d}'.format(self.executor_process.pid))
                os.kill(self.executor_process.pid, signal.SIGKILL)
                self.executor_process.join()
                
            self.log.debug('worker process {:d} terminated with code {:d}'.format(self.executor_process.pid, self.executor_process.exitcode))
        else:
            self.log.debug('worker process {:d} terminated gracefully with code {:d}'.format(self.executor_process.pid, self.executor_process.exitcode))        
        assert not self.executor_process.is_alive()
        

    def startup(self):
        executor = ZMQExecutor(self.task_endpoint, self.result_endpoint)
        self.executor_process = multiprocessing.Process(target = executor.startup)
        self.executor_process.start()
        self.context = zmq.Context()
        self.comm_thread = threading.Thread(target=self.comm_loop)
        self.comm_thread.start()

            
class ZMQExecutor(ZMQCore):
    '''The is the component of the ZMQ WM worker that actually executes tasks.
    This is isolated in a separate process and controlled via ZMQ from 
    the ZMQWorker.'''
    
    def __init__(self, task_endpoint, result_endpoint):
        super(ZMQExecutor,self).__init__()
        
        
        self.task_endpoint = task_endpoint
        self.result_endpoint = result_endpoint
        

    def comm_loop(self):
        
        task_socket = self.context.socket(zmq.PULL)
        result_socket = self.context.socket(zmq.PUSH)
        
        task_socket.bind(self.task_endpoint)
        result_socket.connect(self.result_endpoint)
        
        self.log.info('This is {}'.format(self.node_description))
        
        while True:
            msg = self.recv_message(task_socket)
            
            if msg.message == Message.TASK:
                task = msg.payload
                result = task.execute()
                self.send_message(result_socket, Message.RESULT, result)
            elif msg.message == Message.SHUTDOWN:
                break
            
    def startup(self):
        self.context = zmq.Context()
        self.comm_loop()

