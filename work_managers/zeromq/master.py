'''
Created on May 29, 2015

@author: mzwier
'''

from cPickle import HIGHEST_PROTOCOL
from core import (MSG_MASTER_ACK_ALIVE,MSG_MASTER_Q_ALIVE,
                   MSG_SHUTDOWN,MSG_SHUTDOWN_ACK,
                   MSG_TASK_REQUEST,MSG_TASK_RESULT,
                   MSG_WORKER_ACK_EXISTS,MSG_WORKER_Q_EXISTS)

from core import (MASTER_CRASH_TIMEOUT,WORKER_CRASH_TIMEOUT,DEFAULT_STATUS_POLL)
from core import ZMQCore

import gevent
from zmq import green as zmq

from cPickle import HIGHEST_PROTOCOL

class ZMQMaster(ZMQCore):
    def __init__(self):
        super(ZMQMaster,self).__init__()
        self.rr_endpoint = 'tcp://*:23811'
        self.ann_endpoint = 'tcp://*:23812'
        
        # Amount of time (in s) to wait for contact from any worker before
        # exiting. This lets us watch for the case where the master starts up
        # successfully but the workers never do, or all the workers crash but
        # so much, but can mean a few thousand CPU hours on a supercomputer.
        self.worker_hangcheck_timeout = 30
        self.worker_hangcheck_timer = None
        
        # Amount of time (in s) to announce our presence to workers. This lets
        # workers know that the master exists, so they can exit if the master
        # crashes or hangs. This also tells workers to announce themselves to
        # the master, so that the master knows. This should be small for local
        # work and larger for supercomputer work.
        self.beacon_period = 1
        self.beacon_timer = None
                
        self.worker_info = {}
                
        self.contact_established = False
        self.startup_hangcheck_timeout = 10
            
    def taskreq_handler(self):
        while True:
            message = self.rr_socket.recv_pyobj()
            print('MASTER> received {!r}'.format(message))
            self.rr_socket.send_pyobj('ok')
    
    def shutdown_handler(self, signal=None):
        self.ann_socket.send_pyobj(MSG_SHUTDOWN)
        super(ZMQMaster,self).shutdown_handler(signal)
                
    def comm_loop(self):
        '''Master communication loop for the master process.'''
        
        
        self.context = zmq.Context()
        self.rr_socket = self.context.socket(zmq.REP)
        self.ann_socket = self.context.socket(zmq.PUB)
        
        self.install_signal_handlers()
        taskreq_greenlet = gevent.spawn(self.taskreq_handler)
                 
        try:
            self.rr_socket.bind(self.rr_endpoint)
            self.ann_socket.bind(self.ann_endpoint)
            taskreq_greenlet.join()
        except KeyboardInterrupt:
            self.shutdown_handler()
        finally:
            self.context.destroy(linger=1)
            
        
if __name__ == '__main__':
    zm = ZMQMaster()
    zm.comm_loop()
    