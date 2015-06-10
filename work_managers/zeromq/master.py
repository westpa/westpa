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

class ZMQMaster(ZMQCore):
        
    def __init__(self, upstream_task_endpoint, upstream_result_endpoint, upstream_ann_endpoint,
                 downstream_rr_endpoint, downstream_ann_endpoint):
        super(ZMQMaster,self).__init__()
        
        # Our downstream connections
        # We use request/reply for two reasons:
        #   1.  we can use the same socket to distribute both tasks and results
        #       (with push/pull we'd need one for each direction); this gives us
        #       more nodes before we have to coalesce communications
        #   2.  we anticipate active load-balancing, which requires req/rep to overcome
        #       the default load balancing of push/pull
        self.downstream_rr_endpoint = downstream_rr_endpoint
        self.downstream_ann_endpoint = downstream_ann_endpoint
        self.downstream_rr_socket = None
        self.downstream_ann_socket = None
        
        # Our upstream connection
        # The work manager streams tasks to us with PUSH and retrieves results with PULL
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_result_endpoint = upstream_result_endpoint
        self.upstream_ann_endpoint = upstream_ann_endpoint
        
        # Tasks waiting for dispatch
        self.outgoing_tasks = deque()
        
        # Tasks being processed by workers (indexed by worker_id)
        self.pending_tasks = dict()
                        
        # Worker information, indexed by worker_id
        self.worker_info = {}
        
        # We are the master, so master_id is node_id
        self.master_id = self.node_id
        
        self.comm_thread = None
                
    def send_message(self, socket, message, payload=None, flags=0):
        message = Message(message, payload)
        message.master_id = self.node_id
        super(ZMQMaster,self).send_message(socket, message, payload, flags)
    
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
    
    
    def send_shutdown_message(self, signal=None):
        shutdown_msg = Message(Message.SHUTDOWN, payload=signal, master_id=self.node_id)
        
        # downstream_ann_socket can be None if we are in unit tests
        if self.downstream_ann_socket is not None:
            self.downstream_ann_socket.send_pyobj(shutdown_msg)
        
                
    def comm_loop(self):
        '''Master communication loop for the master process.'''
                
        self.log.info('This is {}'.format(self.node_description))
        
        try:
            self.downstream_rr_socket = self.context.socket(zmq.REP)
            self.downstream_rr_socket.bind(self.downstream_rr_endpoint)
            
            self.downstream_ann_socket = self.context.socket(zmq.PUB)
            self.downstream_ann_socket.bind(self.downstream_ann_endpoint)
            
            self.upstream_ann_socket = self.context.socket(zmq.SUB)
            self.upstream_ann_socket.setsockopt(zmq.SUBSCRIBE, '')
            self.upstream_ann_socket.connect(self.upstream_ann_endpoint)
            
            self.upstream_task_socket = self.context.socket(zmq.PULL)
            self.upstream_task_socket.connect(self.upstream_task_endpoint)
            
            self.upstream_result_socket = self.context.socket(zmq.PUSH)
            self.upstream_result_socket.connect(self.upstream_result_endpoint)
            
            inproc_socket = self.context.socket(zmq.SUB)
            inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
            inproc_socket.bind(self.inproc_endpoint)
            
            poller = zmq.Poller()
            poller.register(inproc_socket, zmq.POLLIN)
            poller.register(self.upstream_ann_socket, zmq.POLLIN)
            poller.register(self.downstream_rr_socket, zmq.POLLIN)
            poller.register(self.upstream_task_socket, zmq.POLLIN)
            
            
            # The communication loop goes a little like this:
            #  1. Check to see if we need to shut down.
            #  2. Accept task requests from upstream and enqueue
            #  3. Respond to worker requests in order.
            #  4. Check timeouts, send heartbeats/tasks avail/etc.
            
            while True:                
                poll_results = dict(poller.poll())
                #log.debug('poll results: {!r}'.format(poll_results))
                
                # Check to see if we need to shut down
                if inproc_socket in poll_results:
                    msgs = self.recv_all(inproc_socket,validate=False)
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        return
                
                if self.upstream_ann_socket in poll_results:
                    msgs = self.recv_all(self.upstream_ann_socket,validate=False)
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        return
                            
                if self.upstream_task_socket in poll_results:
                    msgs = self.recv_all(self.upstream_task_socket)
                    for msg in msgs:
                        log.debug('received {!r}'.format(msg))
                        
                if self.downstream_rr_socket in poll_results:
                    msgs = self.recv_all(self.downstream_rr_socket)
                    for msg in msgs:
                        log.debug('received {!r}'.format(msg))
                        if msg.message == Message.IDENTIFY:
                            self.handle_identify(self.downstream_rr_socket, msg)
                        elif msg.message == Message.TASK_REQUEST:
                            self.handle_task_request(self.downstream_rr_socket, msg)
                        else:
                            self.send_ack(self.downstream_rr_socket, msg)

        finally:
            self.send_shutdown_message()
            self.context.destroy(linger=1)
                        
    def send_tasks_available(self):
        self.send_message(self.downstream_ann_socket,Message.TASKS_AVAILABLE)

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
    