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

import zmq.green as zmq
from cPickle import HIGHEST_PROTOCOL

class ZMQMaster:
    def __init__(self):
        self.rr_endpoint = 'tcp://*:23811'
        self.ann_endpoint = 'tcp://*:23812'
        
        # Amount of time (in s) to wait for contact from any worker before
        # exiting. This lets us watch for the case where the master starts up
        # successfully but the workers never do, or all the workers crash but
        # so much, but can mean a few thousand CPU hours on a supercomputer.
        self.worker_hangcheck_timeout = 30
        
        # Amount of time (in s) to announce our presence to workers. This lets
        # workers know that the master exists, so they can exit if the master
        # crashes or hangs. This also tells workers to announce themselves to
        # the master, so that the master knows
        self.beacon_period = 1
                
        self.worker_info = {}
        
        self.context = None
        self.ann_socket = None
        self.rr_socket = None
        
        
            
    def taskreq_handler(self, socket, message):
        self.contact_established = True
        print('MASTER> received {!r}'.format(message))
        socket.send_pyobj('ok')
    
        
    def comm_loop(self):
        '''Master communication loop for the master process.'''
        
        
        context = zmq.Context()
        rr_socket = context.socket(zmq.REP)
        ann_socket = context.socket(zmq.PUB)
        poller = zmq.Poller()
        
        try:
            rr_socket.bind(self.rr_endpoint)
            ann_socket.bind(self.ann_endpoint)
            
            ann_socket.hwm = 10
            
            poller.register(rr_socket,zmq.POLLIN)
            
            while True:
                poll_results = dict(poller.poll(timeout=(None if self.contact_established else self.startup_hangcheck_timeout * 1000)))
                
                if rr_socket in poll_results:
                    self.taskreq_handler(rr_socket, rr_socket.recv_pyobj())                    
                elif not poll_results:
                    assert not self.contact_established
                    print('MASTER> poll timed out; exiting')
                    break
        except KeyboardInterrupt:
            print('Exiting!\n')
            try:
                ann_socket.send_pyobj(MSG_SHUTDOWN)
            finally:
                pass
            raise
        finally:
            context.destroy(linger=1)
            
if __name__ == '__main__':
    zm = ZMQMaster()
    zm.comm_loop()
    
    