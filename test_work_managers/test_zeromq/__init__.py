import logging
logging.basicConfig(level=logging.DEBUG)

import time, tempfile, os, socket
from unittest import skip

import zmq

from work_managers.zeromq.core import ZMQCore

# Amount of time to wait after executing setUp() to allow sockets to settle
SETUP_WAIT = 0.010

# Amount of time to wait prior to executing tearDown(), to ensure that shutdown
# message don't get lost.
# The original value here (0.010 s = 10 ms) is probably quite generous. 
TEARDOWN_WAIT = 0.010

# How long to wait to let shutdown signals sort themselves out
SHUTDOWN_WAIT = 1

BEACON_PERIOD = 0.2
BEACON_WAIT = BEACON_PERIOD * 5

def sockdelay():
    '''Delay for slightly longer than the default auto-reconnect time for ZeroMQ (100 ms)'''
    time.sleep(0.2)
    
def randport():
    s = socket.socket()
    s.bind(('127.0.0.1',0))
    port = s.getsockname()[1]
    s.close()
    return port

def randipc():
    (fd, socket_path) = tempfile.mkstemp()
    os.close(fd)
    endpoint = 'ipc://{}'.format(socket_path)
    return endpoint

class ZMQTestBase(object):
    '''Support routines'''
    
    # default endpoint type for tests whose transport is not otherwise specified
    endpoint_type = 'ipc' 
    
    def make_ipc_endpoint(self):
        endpoint = randipc()
        try:
            self._endpoints.append(endpoint)
        except AttributeError:
            self._endpoints = [endpoint]
        return endpoint
    
    def make_tcp_endpoint(self):
        return 'tcp://127.0.0.1:{}'.format(randport())
    
    def make_endpoint(self):
        try:
            endpoint_type = self.endpoint_type
        except AttributeError:
            endpoint_type = 'ipc'
        
        if endpoint_type == 'ipc':
            return self.make_ipc_endpoint()
        elif endpoint_type == 'tcp':
            return self.make_tcp_endpoint()
        else:
            raise ValueError('invalid endpoint type set')
            
    def cleanup_endpoints(self):
        for endpoint in self._endpoints:
            try:
                os.unlink(endpoint[6:])
            except OSError:
                pass
        del self._endpoints
        
    def setUp(self):
        self._endpoints = []
        self.test_core = ZMQCore()
        self.test_core.context = self.test_context = zmq.Context()
        self.test_core.validation_fail_action = 'raise'
        
    def tearDown(self):
        self.cleanup_endpoints()
        self.test_context.destroy(linger=1)
        del self.test_context
        
