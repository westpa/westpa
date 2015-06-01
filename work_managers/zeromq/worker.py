'''
Created on May 29, 2015

@author: mzwier
'''

from core import (MSG_MASTER_ACK_ALIVE,MSG_MASTER_Q_ALIVE,
                   MSG_SHUTDOWN,MSG_SHUTDOWN_ACK,
                   MSG_TASK_REQUEST,MSG_TASK_RESULT,
                   MSG_WORKER_ACK_EXISTS,MSG_WORKER_Q_EXISTS)

from core import (MASTER_CRASH_TIMEOUT,WORKER_CRASH_TIMEOUT,DEFAULT_STATUS_POLL)


import zmq.green as zmq

class ZMQWorker:
    def __init__(self):
        self.rr_endpoint = 'tcp://localhost:23811'
        self.ann_endpoint = 'tcp://localhost:23812'
        
            
    def announcement_handler(self, socket, message):
        print('WORKER> received announcement {!r}'.format(message))
        
        if message == MSG_SHUTDOWN:
            print('WORKER> shutting down')
            raise KeyboardInterrupt
        
    def comm_loop(self):
        '''Master communication loop for the worker process.'''
        
        
        context = zmq.Context()
        rr_socket = context.socket(zmq.REQ)
        ann_socket = context.socket(zmq.SUB)
        poller = zmq.Poller()
        
        
        try:
            rr_socket.connect(self.rr_endpoint)
            ann_socket.connect(self.ann_endpoint)
            ann_socket.setsockopt(zmq.SUBSCRIBE,'')
                        
            poller.register(rr_socket,zmq.POLLIN)
            poller.register(ann_socket,zmq.POLLIN)
            
            while True:
                poll_results = dict(poller.poll(timeout=self.status_poll_time*1000))
                
                if ann_socket in poll_results:
                    msg = ann_socket.recv_pyobj()
                    print('WORKER> received announcement {!r}'.format(msg))
                    
                    if msg == MSG_SHUTDOWN:
                        print('WORKER> shutting down')
                        break

        except KeyboardInterrupt:
            print('Exiting!\n')
            raise
        finally:
            context.destroy(linger=1)
            
if __name__ == '__main__':
    zw = ZMQWorker()
    zw.comm_loop()