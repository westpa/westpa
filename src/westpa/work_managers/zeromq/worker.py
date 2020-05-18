'''
Created on May 29, 2015

@author: mzwier
'''

import sys

import logging
from _ast import Break
log = logging.getLogger(__name__)

from .core import ZMQCore, Message, ZMQWMTimeout, PassiveMultiTimer, Task, Result, TIMEOUT_MASTER_BEACON
import threading, multiprocessing, os, signal
from contextlib import contextmanager


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

    @property
    def is_master(self):
        return False

                    
    def update_master_info(self, msg):
        if self.master_id is None:
            self.master_id = msg.master_id
        self.timers.reset(TIMEOUT_MASTER_BEACON)
        
    def identify(self, rr_socket):
        if self.master_id is None or self.identified or self.timers.expired(TIMEOUT_MASTER_BEACON): return
        self.send_message(rr_socket, Message.IDENTIFY, payload=self.get_identification())
        self.recv_ack(rr_socket,timeout=self.master_beacon_period*self.timeout_factor*1000)
        self.identified = True
        
    def request_task(self, rr_socket, task_socket):
        if self.master_id is None: return
        elif self.pending_task is not None: return
        elif self.timers.expired(TIMEOUT_MASTER_BEACON): return
        else:
            self.send_message(rr_socket, Message.TASK_REQUEST)
            reply = self.recv_message(rr_socket,timeout=self.master_beacon_period*self.timeout_factor*1000)
            self.update_master_info(reply)
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
        reply = self.recv_ack(rr_socket, timeout=self.master_beacon_period*self.timeout_factor*1000)
        self.update_master_info(reply)
            
    def comm_loop(self):
        '''Master communication loop for the worker process.'''
        
        rr_socket = self.context.socket(zmq.REQ)

        ann_socket = self.context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,b'')
        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,b'')
        
        task_socket = self.context.socket(zmq.PUSH)
        result_socket = self.context.socket(zmq.PULL)
        
        self.log.info('This is {}'.format(self.node_description))
        
        timers = self.timers = PassiveMultiTimer()
        timers.add_timer(TIMEOUT_MASTER_BEACON, 86400)
        timers.add_timer('worker_beacon', self.worker_beacon_period)
        timers.add_timer('startup_timeout', self.startup_timeout)
        peer_found = False
        
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
                
                if poll_results and not peer_found:
                    timers.remove_timer('startup_timeout')
                    peer_found = True
                    timers.change_duration(TIMEOUT_MASTER_BEACON, self.master_beacon_period*self.timeout_factor)
                    timers.reset(TIMEOUT_MASTER_BEACON)
                
                announcements = []   
                
                # Check for internal messages first
                if inproc_socket in poll_results:
                    announcements.extend(self.recv_all(inproc_socket))
          
                # Process announcements
                if ann_socket in poll_results:
                    announcements.extend(self.recv_all(ann_socket))
                    
                #announcements = Message.coalesce_announcements(announcements)
                #self.log.debug('received {:d} announcements'.format(len(announcements)))
                
                messages_by_tag = {}
                for msg in announcements:
                    messages_by_tag.setdefault(msg.message, list()).append(msg)
                    
                # Check for shutdown messages
                if Message.SHUTDOWN in messages_by_tag:
                    self.log.debug('received shutdown message')
                    return
                elif Message.TASKS_AVAILABLE in messages_by_tag:
                    self.update_master_info(messages_by_tag[Message.TASKS_AVAILABLE][0])
                elif Message.MASTER_BEACON in messages_by_tag:
                    self.update_master_info(messages_by_tag[Message.MASTER_BEACON][0])
                    
                if self.master_id is not None and timers.expired('worker_beacon'):
                    self.identify(rr_socket)
                    timers.reset('worker_beacon')
                
                # Handle results, so that we clear ourselves of completed tasks
                # before asking for more
                if result_socket in poll_results:
                    self.handle_result(result_socket, rr_socket)
                    # immediately request another task if available
                    if not timers.expired(TIMEOUT_MASTER_BEACON):
                        self.request_task(rr_socket, task_socket)
                
                # Handle any remaining messages
                for tag, msgs in messages_by_tag.items():
                    if tag == Message.MASTER_BEACON:
                        self.identify(rr_socket)
                    elif tag == Message.RECONFIGURE_TIMEOUT:
                        for msg in msgs:
                            self.handle_reconfigure_timeout(msg, timers)
                    elif tag == Message.TASKS_AVAILABLE:
                        self.request_task(rr_socket,task_socket)      
                        
                del announcements, messages_by_tag
            
                if timers.expired(TIMEOUT_MASTER_BEACON):
                    self.log.error('no contact from master; shutting down')
                    return
                
                if not peer_found and timers.expired('startup_timeout'):
                    self.log.error('startup phase elapsed with no contact from master; shutting down')
                    break
        
        except ZMQWMTimeout:
            # both handle_result() and request_task() have receive timeouts set to 
            # self.master_beacon_period*self.timeout_factor, and timeout exceptions
            # propagate to this point.
            self.log.error('timeout communicating with peer; shutting down')
        finally:
            self.shutdown_executor()
            self.executor_process.join()
            self.context.destroy(linger=1)
            self.context = None
            self.remove_ipc_endpoints()
            
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

        pid = self.executor_process.pid
        self.executor_process.join(self.shutdown_timeout)
        # is_alive() is prone to a race condition so catch the case that the PID is already dead
        if self.executor_process.is_alive():
            self.log.debug('sending SIGTERM to worker process {:d}'.format(pid))
            self.executor_process.terminate()
            # try:
            #     os.kill(self.executor_process.pid, signal.SIGINT)
            # except ProcessLookupError:
            #     self.log.debug('worker process {:d} already dead'.format(pid))
            self.executor_process.join(self.shutdown_timeout)
            if self.executor_process.is_alive():
                self.executor_process.kill()
                self.log.warning('sending SIGKILL to worker process {:d}'.format(pid))
                # try:
                #     os.kill(self.executor_process.pid, signal.SIGKILL)
                # except ProcessLookupError:
                #     self.log.debug('worker process {:d} already dead'.format(pid))
            self.executor_process.join()
            self.log.debug('worker process {:d} terminated'.format(pid))
        else:
            self.log.debug('worker process {:d} terminated gracefully with code {:d}'.format(pid, self.executor_process.exitcode))
        
    def install_signal_handlers(self, signals = None):
        if not signals:
            signals = {signal.SIGINT, signal.SIGQUIT, signal.SIGTERM}
        
        for sig in signals:
            signal.signal(sig, signal.SIG_IGN)

    def startup(self, process_index=None):
        self.install_signal_handlers()
        executor = ZMQExecutor(self.task_endpoint, self.result_endpoint)
        self.executor_process = multiprocessing.Process(target = executor.startup, args=(process_index,))
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
        
        try:
            while True:
                try:
                    msg = self.recv_message(task_socket,timeout=100)
                except KeyboardInterrupt:
                    break
                except ZMQWMTimeout:
                    continue
                else:
                    if msg.message == Message.TASK:
                        task = msg.payload
                        result = task.execute()
                        self.send_message(result_socket, Message.RESULT, result)
                    elif msg.message == Message.SHUTDOWN:
                        break
        finally:
            self.context.destroy(linger=0)
            self.context = None
            
            
    def startup(self, process_index=None):
        if process_index is not None:
            from work_managers import environment
            pi_name = '{}_PROCESS_INDEX'.format(environment.WMEnvironment.env_prefix)
            self.log.debug('Setting {}={}'.format(pi_name, process_index))
            os.environ[pi_name] = str(process_index)
        self.context = zmq.Context()
        self.comm_loop()

