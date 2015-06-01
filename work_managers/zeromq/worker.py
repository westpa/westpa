'''
Created on May 29, 2015

@author: mzwier
'''

from core import (MSG_MASTER_ACK_ALIVE,MSG_MASTER_Q_ALIVE,
                   MSG_SHUTDOWN,MSG_SHUTDOWN_ACK,
                   MSG_TASK_REQUEST,MSG_TASK_RESULT,
                   MSG_WORKER_ACK_EXISTS,MSG_WORKER_Q_EXISTS)

from core import (MASTER_CRASH_TIMEOUT,WORKER_CRASH_TIMEOUT,DEFAULT_STATUS_POLL)
from core import ZMQCore

import gevent
import zmq.green as zmq

class ZMQWorker(ZMQCore):
    def __init__(self):
        super(ZMQWorker,self).__init__()
        self.rr_endpoint = 'tcp://localhost:23811'
        self.ann_endpoint = 'tcp://localhost:23812'
        
        self.status_poll_time = 2
        
            
    def announcement_handler(self):
        while True:
            message = self.ann_socket.recv_pyobj()
            print('WORKER> received {!r}'.format(message))
            if message == MSG_SHUTDOWN:
                print('WORKER> shutting down')
                return
            
    def taskreq_handler(self):
        for n in xrange(5):
            self.rr_socket.send_pyobj('hello?')
            print('WORKER> {!r}'.format(self.rr_socket.recv_pyobj()))
            gevent.sleep(2)
        
    def comm_loop(self):
        '''Master communication loop for the worker process.'''
        
        
        self.context = zmq.Context()
        self.rr_socket = self.context.socket(zmq.REQ)
        self.ann_socket = self.context.socket(zmq.SUB)
        
        
        try:
            self.rr_socket.connect(self.rr_endpoint)
            self.ann_socket.connect(self.ann_endpoint)
            self.ann_socket.setsockopt(zmq.SUBSCRIBE,'')
            
            task_greenlet = gevent.spawn(self.taskreq_handler)
            ann_greenlet = gevent.spawn(self.announcement_handler)
            ann_greenlet.join()
            
        finally:
            self.context.destroy(linger=1)
            
if __name__ == '__main__':
    zw = ZMQWorker()
    zw.comm_loop()