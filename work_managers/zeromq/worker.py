'''
Created on May 29, 2015

@author: mzwier
'''

import sys

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, ZMQWMEnvironmentError

#import zmq
import zmq.green as zmq

class ZMQWorker(ZMQCore):
    '''This is the outward facing worker component of the ZMQ work manager. This
    forms the interface to the master. This process cannot hang or crash due to an 
    error in tasks it executes, so tasks are isolated in ZMQExecutor, which 
    communicates with ZMQWorker via (what else?) ZeroMQ.'''
    
    def __init__(self, rr_endpoint, ann_endpoint):
        super(ZMQWorker,self).__init__()
        self.rr_endpoint = rr_endpoint
        self.ann_endpoint = ann_endpoint
        self.master_id = None
        self.last_master_contact = None
                
    def validate_message(self, message):
        '''Validate an incoming message'''
        super(ZMQWorker,self).validate_message(message)
        if self.master_id is not None and message.master_id != self.master_id:
            raise ZMQWMEnvironmentError('incoming message expected from master {!s} but obtained from {!s}'
                                        .format(self.master_id, message.master_id))
        if message.dest_id is not None and message.dest_id != self.node_id:
            raise ZMQWMEnvironmentError('incoming message destined for node {!s} but this is node {!s}'
                                        .format(message.dest_id, self.node_id))                      


    def send_message(self, socket, message, payload=None, flags=0):
        message = Message(message, payload)
        message.master_id = self.master_id
        message.worker_id = self.node_id
        super(ZMQWorker,self).send_message(socket, message,flags)
        
    def handle_pairing(self):
        self.send_message(self.rr_socket, Message.IDENTIFY, payload=self.get_identification())
        initial_reply = self.recv_message(self.rr_socket)
        if self.master_id is None:
            self.master_id = initial_reply.master_id
        self.log.debug('server identity: {!r}'.format(initial_reply.payload))
        self.log.info('paired with master {!s}'.format(self.master_id))
        
    def comm_loop(self):
        '''Master communication loop for the worker process.'''
        
        self.rr_socket = self.context.socket(zmq.REQ)

        self.ann_socket = self.context.socket(zmq.SUB)
        inproc_socket = self.context.socket(zmq.SUB)
        
        self.log.info('This is {}'.format(self.node_description))
        
       
        
        try:
            self.rr_socket.connect(self.rr_endpoint)
            
            self.ann_socket.connect(self.ann_endpoint)
            self.ann_socket.setsockopt(zmq.SUBSCRIBE,'')
            

            inproc_socket.bind(self.inproc_endpoint)            
            inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
            
            poller = zmq.Poller()
            poller.register(self.ann_socket, zmq.POLLIN)
            poller.register(inproc_socket, zmq.POLLIN)
            
            while True: 
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
                                    
          
                    
                if self.ann_socket in poll_results:
                    while True:
                        try:
                            msg = self.recv_message(self.ann_socket,zmq.NOBLOCK)
                        except zmq.Again:
                            break
                        else:
                            self.log.info('received {!r}'.format(msg))
                            if msg.message == Message.TASKS_AVAILABLE:
                                self.log.debug('tasks available')
                            elif msg.message == Message.SHUTDOWN:
                                self.log.info('shutting down')
                                return
                            elif msg.message == Message.MASTER_BEACON:
                                pass
                        
                            
        finally:
            self.context.destroy(linger=0)
            
class ZMQExecutor(ZMQCore):
    '''The is the component of the ZMQ WM worker that actually executes tasks.
    This is isolated in a separate process and controlled via ZMQ from 
    the ZMQWorker.'''
    
    pass

            
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    
    mode = 'ipc'
    
    if mode == 'ipc':
        rr_endpoint = 'ipc:///tmp/zmqwmtest.rr'
        ann_endpoint = 'ipc:///tmp/zmqwmtest.ann'
    else: # mode == 'tcp':
        rr_endpoint='tcp://localhost:23811'
        ann_endpoint='tcp://localhost:23812'
    
    try:    
        zw = ZMQWorker(rr_endpoint, ann_endpoint)
        zw.validation_fail_action = 'warn'
        zw.startup()
        zw.join()
    finally:
        ZMQWorker.remove_ipc_endpoints()
