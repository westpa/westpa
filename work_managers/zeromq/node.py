'''
Created on Jun 11, 2015

@author: mzwier
'''


import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message

import zmq
from zmq.devices import ThreadProxy

class ZMQNode(ZMQCore):
    def __init__(self, upstream_rr_endpoint, upstream_ann_endpoint):
        super(ZMQNode,self).__init__()
        
        self.upstream_rr_endpoint = upstream_rr_endpoint
        self.upstream_ann_endpoint = upstream_ann_endpoint
        
        self.downstream_rr_endpoints = []
        self.downstream_ann_endpoints = []
        
    
    def comm_loop(self):
        self.context = zmq.Context()
        
        rr_proxy = ThreadProxy(zmq.ROUTER, zmq.DEALER)
        ann_proxy = ThreadProxy(zmq.SUB, zmq.PUB, zmq.PUB)
        
        ann_monitor = self.context.socket(zmq.SUB)
        ann_monitor.setsockopt(zmq.SUBSCRIBE, '')
        
        for endpoint in self.downstream_rr_endpoints:
            rr_proxy.bind_in(endpoint)
        rr_proxy.connect_out(self.upstream_rr_endpoint)
            
        for endpoint in self.downstream_ann_endpoints:
            ann_proxy.bind_out(endpoint)
        ann_proxy.connect_in(self.upstream_ann_endpoint)
        ann_proxy.setsockopt_in(zmq.SUBSCRIBE, '')
        
        #mon_endpoint = 'inproc://fo358728735'
        mon_endpoint = self.make_internal_endpoint()
        ann_proxy.bind_mon(mon_endpoint)
        ann_monitor.connect(mon_endpoint)
        
        rr_proxy.start()
        ann_proxy.start()
        
        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
        inproc_socket.bind(self.inproc_endpoint)        
        
        poller = zmq.Poller()
        poller.register(ann_monitor, zmq.POLLIN)
        poller.register(inproc_socket, zmq.POLLIN)
        try:
            while True:
                poll_results = dict(poller.poll())
                
                if inproc_socket in poll_results:
                    msgs = self.recv_all(ann_monitor,validate=False)
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        self.log.debug('shutdown received')
                        break                    
                
                if ann_monitor in poll_results:
                    msgs = self.recv_all(ann_monitor,validate=False)
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        self.log.debug('shutdown received')
                        break                    
            
        finally:
            self.log.debug('exiting')
            self.context.destroy(linger=100)
            self.context = None
            self.remove_ipc_endpoints()
