'''
Created on Jun 11, 2015

@author: mzwier
'''


import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, PassiveMultiTimer, IsNode

import zmq
from zmq.devices import ThreadProxy

class ZMQNode(ZMQCore,IsNode):
    def __init__(self, upstream_rr_endpoint, upstream_ann_endpoint, n_local_workers=None):
        ZMQCore.__init__(self)
        IsNode.__init__(self,n_local_workers)
        
        self.upstream_rr_endpoint = upstream_rr_endpoint
        self.upstream_ann_endpoint = upstream_ann_endpoint

    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_traceback):
        return False
    
    def run(self):
        self.startup()
            
    @property
    def is_master(self):
        return False
        
        
    def comm_loop(self):
        self.context = zmq.Context.instance() 
        # or else the proxies create sockets in a different context
         
        self.context.linger = 100
        # So we don't have to destroy the context at the end of the loop
        
        rr_proxy = ThreadProxy(zmq.ROUTER, zmq.DEALER)
        
        # We use push/pull so (1) we don't miss any announcements
        # and (2) we don't have to deal with subscription messages
        ann_proxy = ThreadProxy(zmq.SUB, zmq.PUB, zmq.PUSH)
        ann_monitor = self.context.socket(zmq.PULL)
        
        # Not monitoring request/reply streams for two reasons:
        # (1) we'd need to strip identity frames to interpret the messages
        # (2) interpreting the messages means we'd have to decode (unpickle) and then re-encode
        # all of the data flying through here, which seems like a waste just to see if
        # clients start up. We miss the edge failure case where one node's workers
        # start up but another's fail. Seems much less likely than all workers
        # failing to start up, which would be caught by the master
                
        ann_mon_endpoint = 'inproc://{:x}'.format(id(ann_monitor))
        ann_monitor.bind(ann_mon_endpoint)
        
    
        rr_proxy.bind_in(self.downstream_rr_endpoint)
        if self.local_rr_endpoint: rr_proxy.bind_in(self.local_rr_endpoint)
        self.log.debug('connecting upstream_rr_endpoint = {!r}'.format(self.upstream_rr_endpoint))
        rr_proxy.connect_out(self.upstream_rr_endpoint)
            
        ann_proxy.bind_out(self.downstream_ann_endpoint)
        if self.local_ann_endpoint: ann_proxy.bind_out(self.local_ann_endpoint)
        ann_proxy.connect_in(self.upstream_ann_endpoint)
        self.log.debug('connecting upstream_ann_endpoint = {!r}'.format(self.upstream_ann_endpoint))        
        ann_proxy.setsockopt_in(zmq.SUBSCRIBE, '')        
        ann_proxy.connect_mon(ann_mon_endpoint)
        
        rr_proxy.start()
        ann_proxy.start()
        
        ann_monitor.connect(ann_mon_endpoint)
        
        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
        inproc_socket.bind(self.inproc_endpoint)
        
        timers = PassiveMultiTimer()
        timers.add_timer('master_beacon', self.master_beacon_period)
        timers.add_timer('startup_timeout', self.startup_timeout)
        timers.reset()
        
        self.log.debug('master beacon period: {!r}'.format(self.master_beacon_period))
        self.log.debug('startup timeout: {!r}'.format(self.startup_timeout))
        
        peer_found = False
        
        poller = zmq.Poller()
        poller.register(ann_monitor, zmq.POLLIN)
        poller.register(inproc_socket, zmq.POLLIN)
        try:
            while True:
                poll_results = dict(poller.poll((timers.next_expiration_in() or 0.001)*1000))
                
                if inproc_socket in poll_results:
                    msgs = self.recv_all(ann_monitor,validate=False)
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        self.log.debug('shutdown received')
                        break                    
                
                if ann_monitor in poll_results:
                    msgs = self.recv_all(ann_monitor,validate=False)
                    message_tags = {msg.message for msg in msgs}
                    if Message.SHUTDOWN in message_tags:
                        self.log.debug('shutdown received')
                        break
                    if not peer_found and (Message.MASTER_BEACON in message_tags or Message.TASKS_AVAILABLE in message_tags):
                        peer_found = True
                        timers.remove_timer('startup_timeout')
                        
                        
                if not peer_found and timers.expired('startup_timeout'):
                    self.log.error('startup phase elapsed with no contact from peer; shutting down')
                    break
            
        finally:
            self.log.debug('exiting')
            self.context = None
            self.remove_ipc_endpoints()
            IsNode.shutdown(self)

    def startup(self):
        IsNode.startup(self)
        super(ZMQNode,self).startup()
        