'''
Created on May 29, 2015

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, Task, ZMQWMEnvironmentError, ZMQWorkerMissing

#import gevent
#from zmq import green as zmq
import zmq
import time
from collections import deque
import threading

from cPickle import HIGHEST_PROTOCOL

class ZMQMaster(ZMQCore):
        
    def __init__(self, rr_endpoint, ann_endpoint):
        super(ZMQMaster,self).__init__()
        
        # Our downstream connections
        self.rr_endpoint = rr_endpoint
        self.ann_endpoint = ann_endpoint
        self.rr_socket = None
        self.ann_socket = None
        
        # Our upstream connection
        
        # Tasks waiting for dispatch
        self.outgoing_tasks = deque()
        
        # Tasks being processed by workers (indexed by worker_id)
        self.pending_tasks = dict()
                        
        # Worker information, indexed by worker_id
        self.worker_info = {}
        
        # We are the master, so master_id is node_id
        self.master_id = self.node_id
        
        self.comm_thread = None
        
    def validate_message(self, message):
        '''Validate an incoming message'''
        super(ZMQMaster,self).validate_message(message)
        if message.master_id not in (self.node_id, None):
            raise ZMQWMEnvironmentError('incoming message destined for another master (this={!s}, incoming={!s}'.format(self.node_id, message.master_id))
        elif message.worker_id is None:
            raise ZMQWMEnvironmentError('incoming message has blank worker ID')
        
    def send_message(self, socket, message, payload=None):
        message = Message(message, payload)
        message.master_id = self.node_id
        super(ZMQMaster,self).send_message(socket, message)
    
    def update_worker_lastseen(self, msg, seen_at=None):
        worker_record = self.worker_info.setdefault(msg.worker_id, {})
        worker_record['last_seen'] = seen_at or time.time()
        
    def update_worker_identification(self, worker_id, identification):
        worker_record = self.worker_info.setdefault(worker_id, {})
        worker_record['identification'] = identification
        self.log.debug('worker identification updated: {!r}'.format(worker_record))
        
    def drop_worker(self, worker_id):
        self.log.warning('dropping worker {}'.format(worker_id))
        self.worker_info.pop(worker_id, None)
        slayed_task = self.pending_tasks.pop(worker_id, None)
        if slayed_task:
            # Currently, just report an error; in the future, we could re-queue, but 
            # only if tasks are structured so that partial results are overwritten
            slayed_task.future._set_exception(ZMQWorkerMissing('lost contact with worker processing this task'))
        
    def handle_identify(self, socket, msg):
        assert msg.message == Message.IDENTIFY
        
        # Record worker identification
        identification = msg.payload
        with self.message_validation(msg):
            assert isinstance(identification, dict)
            self.update_worker_lastseen(msg)
            self.update_worker_identification(msg.worker_id, identification)
            
        # Reply with our identification
        self.send_reply(socket, msg, Message.IDENTIFY, payload=self.get_identification())
    
    
    def rr_handler(self):
        while True:
            message = self.recv_message(self.rr_socket)
            
    
    def send_shutdown_message(self, signal=None):
        shutdown_msg = Message(Message.SHUTDOWN, payload=signal, master_id=self.node_id)
        self.ann_socket.send_pyobj(shutdown_msg)

    #def shutdown_handler(self, signal=None, frame=None):
    #    self.send_shutdown_message(signal)
    #    super(ZMQMaster,self).shutdown_handler(signal, frame)
        
    def beacon_handler(self):
        while True:
            self.ann_socket.send_pyobj(Message(master_id=self.node_id, message=Message.MASTER_BEACON))
            if self.outgoing_tasks:
                self.send_tasks_available()
            gevent.sleep(self.master_beacon_period)
            
    def worker_error_checker(self):
        while True:
            # This loops no *more often* than once every ``worker_beacon_period`` seconds
            gevent.sleep(self.worker_beacon_period)
            
            now = time.time()
            for worker_id, worker_info in self.worker_info.items():
                if now - worker_info.get('last_seen', now) > self.worker_beacon_period:
                    self.log.warning('worker {} is not responding'.format(worker_id))
                    self.drop_worker(worker_id)
                
    def comm_loop(self):
        '''Master communication loop for the master process.'''
        
        self.rr_socket = self.context.socket(zmq.REP)
        self.ann_socket = self.context.socket(zmq.PUB)        
        inproc_socket = self.context.socket(zmq.SUB)
        
        self.log.info('This is {}'.format(self.node_description))
        
        try:
            inproc_socket.bind(self.inproc_endpoint)
            inproc_socket.setsockopt(zmq.SUBSCRIBE,'')

            self.rr_socket.bind(self.rr_endpoint)
            self.ann_socket.bind(self.ann_endpoint)
            
            poller = zmq.Poller()
            poller.register(self.rr_socket, zmq.POLLIN)
            poller.register(inproc_socket, zmq.POLLIN)
            while True:
                if self.outgoing_tasks:
                    self.send_tasks_available()
                
                poll_results = dict(poller.poll())
                
                # Check for internal messages first
                if inproc_socket in poll_results:
                    while True:
                        try:
                            msg = self.recv_message(inproc_socket,zmq.NOBLOCK,validate=False)
                        except zmq.Again:
                            break
                        else:
                            if msg.message == Message.SHUTDOWN:
                                return
                
                # Read all available messages
                if self.rr_socket in poll_results:
                    while True:
                        try:
                            msg = self.recv_message(self.rr_socket, zmq.NOBLOCK)
                        except zmq.Again:
                            # No more messages to process at this time
                            break
                        else:
                            self.update_worker_lastseen(msg)
                            self.log.info('received {!r}'.format(msg))
                            
                            if msg.message == Message.IDENTIFY:
                                self.handle_identify(self.rr_socket, msg)
                            elif msg.message == Message.TASK_REQUEST:
                                self.handle_task_request(self.rr_socket, msg)
                            else:
                                self.send_ack(self.rr_socket, msg)                
        finally:
            self.send_shutdown_message()
            self.context.destroy(linger=1)
                        
    def send_tasks_available(self):
        self.send_message(self.ann_socket,Message.TASKS_AVAILABLE)

#     def submit(self, fn, args=None, kwargs=None):
#         return self.submit_many([(fn,
#                                   args if args is not None else (),
#                                   kwargs if kwargs is not None else {})]
#                                 )[0]
#     
#     def submit_many(self, tasks):
#         futures = []
#         outgoing_tasks = self.outgoing_tasks
#         
#         for (fn,args,kwargs) in tasks:
#             task = Task(fn, args, kwargs)
#             outgoing_tasks.append(task)
#             futures.append(task.future)
#                     
#         return futures
            
        
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    
    from test_work_managers.tsupport import identity
    
    mode = 'ipc'
    
    if mode == 'ipc':
        rr_endpoint = 'ipc:///tmp/zmqwmtest.rr'
        ann_endpoint = 'ipc:///tmp/zmqwmtest.ann'
        ZMQMaster._ipc_endpoints_to_delete.extend([rr_endpoint,ann_endpoint])
    else: # mode == 'tcp':
        rr_endpoint='tcp://*:23811'
        ann_endpoint='tcp://*:23812'
    
    try:    
        zm = ZMQMaster(rr_endpoint, ann_endpoint)
        zm.validation_fail_action = 'warn'
        #zm.comm_loop()
        #zm.install_signal_handlers()
        zm.startup()
        zm.submit(identity,(1,))
        zm.join()
    finally:
        ZMQMaster.remove_ipc_endpoints()
    