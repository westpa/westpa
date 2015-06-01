'''
Created on May 29, 2015

@author: mzwier
'''

from core import ZMQCore, Message, ZMQWMEnvironmentError

import gevent
import zmq.green as zmq

class ZMQWorker(ZMQCore):
    def __init__(self):
        super(ZMQWorker,self).__init__()
        self.rr_endpoint = 'tcp://localhost:23811'
        self.ann_endpoint = 'tcp://localhost:23812'
        self.master_id = None
        
    def send_message(self, socket, message=None, payload=None):
        socket.send_pyobj(Message(message, payload, master_id = self.master_id, worker_id = self.id))
        
    def validate_message(self, message):
        '''Validate an incoming message'''
        super(ZMQWorker,self).validate_message(message)
        if self.master_id is not None and message.master_id != self.master_id:
            raise ZMQWMEnvironmentError('incoming message expected from master {!s} but obtained from {!s}'
                                        .format(self.master_id, message.master_id))
        if message.worker_id is not None and message.worker_id != self.id:
            raise ZMQWMEnvironmentError('incoming message destined for worker {!s} but this is worker {!s}'
                                        .format(message.worker_id, self.id))                      
            
    def announcement_handler(self):
        while True:
            message = self.recv_message(self.ann_socket)
            self.dumpinfo('received {!r}'.format(message))
            if message.message == Message.SHUTDOWN:
                self.dumpinfo('shutting down')
                return
            elif message.message == Message.MASTER_BEACON:
                pass
                    
    def taskreq_handler(self):
        # First, identify server. Do nothing until we have.
        
        self.rr_socket.send_pyobj(Message(Message.IDENTIFY_MASTER, worker_id=self.id))
        initial_reply = self.recv_message(self.rr_socket)
        if self.master_id is None:
            self.master_id = initial_reply.master_id
        self.dumpinfo('paired with master {!s}'.format(self.master_id))    
        
        for _n in xrange(5):
            self.send_message(self.rr_socket, 'hello?')
            self.recv_message(self.rr_socket)
            gevent.sleep(2)
        
    def comm_loop(self):
        '''Master communication loop for the worker process.'''
        
        self.context = zmq.Context()
        self.rr_socket = self.context.socket(zmq.REQ)
        self.ann_socket = self.context.socket(zmq.SUB)
        
        self.install_signal_handlers()
        
        try:
            self.rr_socket.connect(self.rr_endpoint)
            self.ann_socket.connect(self.ann_endpoint)
            
            # We exit upon receiving a shutdown announcement, so we join the greenlet that's
            # listening for announcements, which exits when a shutdown message is received.
            self.ann_socket.setsockopt(zmq.SUBSCRIBE,'')
            
            task_greenlet = gevent.spawn(self.taskreq_handler)
            ann_greenlet = gevent.spawn(self.announcement_handler)
            ann_greenlet.join()
            
        finally:
            self.context.destroy(linger=1)
            
if __name__ == '__main__':
    zw = ZMQWorker()
    zw.comm_loop()