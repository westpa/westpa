'''
Created on May 29, 2015

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, ZMQWMEnvironmentError

import gevent
from zmq import green as zmq
import time

from cPickle import HIGHEST_PROTOCOL

class ZMQMaster(ZMQCore):
        
    def __init__(self, rr_endpoint, ann_endpoint):
        super(ZMQMaster,self).__init__()
        self.rr_endpoint = rr_endpoint
        self.ann_endpoint = ann_endpoint
        
        # Amount of time (in s) to wait for contact from any worker before
        # exiting. This lets us watch for the case where the master starts up
        # successfully but the workers never do, or all the workers crash but
        # the master doesn't. This doesn't mean too much locally, but can save a
        # few thousand CPU hours on a supercomputer.
        self.worker_hangcheck_timeout = 10
        
        # Amount of time (in s) to announce our presence to workers. This lets
        # workers know that the master exists, so they can exit if the master
        # crashes or hangs. This also tells workers to announce themselves to
        # the master, so that the master knows. This should be small for local
        # work and larger for supercomputer work.
        self.beacon_period = 5
                
        self.worker_info = {}
                
        self.master_id = self.node_id
        
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
        # drop tasks too
        self.log.warning('dropping worker {}'.format(worker_id))
        self.worker_info.pop(worker_id, None)
        
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
        
    def taskreq_handler(self):
        while True:
            message = self.recv_message(self.rr_socket)
            self.update_worker_lastseen(message)
            self.log.info('received {!r}'.format(message))
            
            if message.message == Message.IDENTIFY:
                self.handle_identify(self.rr_socket, message)
            else:
                self.rr_socket.send_pyobj(Message(Message.ACK, master_id=self.node_id, worker_id=message.worker_id))
    
    def send_shutdown_message(self, signal=None):
        shutdown_msg = Message(Message.SHUTDOWN, payload=signal, master_id=self.node_id)
        self.ann_socket.send_pyobj(shutdown_msg)

    def shutdown_handler(self, signal=None):
        self.send_shutdown_message(signal)
        super(ZMQMaster,self).shutdown_handler(signal)
        
    def beacon_handler(self):
        while True:
            self.ann_socket.send_pyobj(Message(master_id=self.node_id, message=Message.MASTER_BEACON))
            gevent.sleep(self.beacon_period)
            
    def worker_error_checker(self):
        while True:
            # This loops no *more often* than once every ``worker_hangcheck_timeout`` seconds
            gevent.sleep(self.worker_hangcheck_timeout)
            
            now = time.time()
            for worker_id, worker_info in self.worker_info.items():
                if now - worker_info.get('last_seen', now) > self.worker_hangcheck_timeout:
                    self.log.warning('worker {} is not responding'.format(worker_id))
                    self.drop_worker(worker_id)
                
    def comm_loop(self):
        '''Master communication loop for the master process.'''
        
        
        self.context = zmq.Context()
        self.rr_socket = self.context.socket(zmq.REP)
        self.ann_socket = self.context.socket(zmq.PUB)
        
        self.install_signal_handlers()
        
        self.log.info('This is {}'.format(self.node_description))
        
        try:
            self.rr_socket.bind(self.rr_endpoint)
            self.ann_socket.bind(self.ann_endpoint)
            taskreq_greenlet = gevent.spawn(self.taskreq_handler)
            beacon_greenlet = gevent.spawn(self.beacon_handler)
            hangcheck_greenlet = gevent.spawn(self.worker_error_checker)
            taskreq_greenlet.join()
        finally:
            self.send_shutdown_message()
            self.context.destroy(linger=1)
            
        
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    zm = ZMQMaster(rr_endpoint='tcp://*:23811',ann_endpoint='tcp://*:23812')
    zm.validation_fail_action = 'warn'
    zm.comm_loop()
    