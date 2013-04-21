# Copyright (C) 2013 Matthew C. Zwier
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

'''
Common constants, objects, routines, and base classes for the ZeroMQ work
manager.
'''

from __future__ import division, print_function; __metaclass__ = type

import os, socket, threading, uuid, tempfile, errno, logging
import zmq

log = logging.getLogger(__name__)

from . import randport, DEFAULT_SERVER_HEARTBEAT_INTERVAL

# General message format: tuples of the form (tag, server_id, client_id, payload),
# where tag is a string, server_id is a UUID or None, client_id is a UUID or 
# None, and payload is any associated data

###################
# Control messages
MSG_ACK = 'ok'
MSG_SHUTDOWN = 'shutdown'
MSG_PING = 'ping' # ping inquiry
MSG_PONG = 'pong' # ping reply
MSG_TASK_REQUEST = 'task'
MSG_TASK_AVAILABLE = 'task_avail'
MSG_TASK_UNAVAILABLE = 'task_unavail'
MSG_RESULT_SUBMISSION = 'result'

# Result types
RESULT_TYPE_RETVAL = 'retval'
RESULT_TYPE_EXCEPTION = 'exception'


class ZMQBase:
    _ipc_endpoints = []

    def __init__(self, server_heartbeat_interval = DEFAULT_SERVER_HEARTBEAT_INTERVAL):
        # ZeroMQ context
        self.context = None
        
        # number of seconds between announcements of where to connect to the master
        self.server_heartbeat_interval = server_heartbeat_interval
        
        # This hostname
        self.hostname = socket.gethostname()
        self.host_id = '{:s}-{:d}'.format(self.hostname, os.getpid())
        self.instance_id = uuid.uuid4()
        
        self._tls = threading.local() 

    @staticmethod
    def canonicalize_endpoint(endpoint, allow_wildcard_host = True):
        if endpoint.startswith('ipc://'):
            return endpoint
        elif endpoint.startswith('tcp://'):
            fields = endpoint[6:].split(':')
            
            # get IP address
            if fields[0] != '*':
                ipaddr = socket.gethostbyname(fields[0])
            else:
                if allow_wildcard_host:
                    ipaddr = '*'
                else:
                    raise ValueError('wildcard host not permitted')
            
            # get/generate port
            try:
                port = fields[1]
            except IndexError:
                # no port given; select one
                port = randport()
            else:
                port = int(fields[1])
                
            return 'tcp://{}:{}'.format(ipaddr,port)
        else:
            raise ValueError('unrecognized/unsupported endpoint: {!r}'.format(endpoint))

    @classmethod    
    def make_ipc_endpoint(cls):
        (fd, socket_path) = tempfile.mkstemp()
        os.close(fd)
        endpoint = 'ipc://{}'.format(socket_path)
        cls._ipc_endpoints.append(endpoint)
        return endpoint
    
    @classmethod
    def remove_ipc_endpoints(cls):
        while cls._ipc_endpoints:
            endpoint = cls._ipc_endpoints.pop()
            assert endpoint.startswith('ipc://')
            socket_path = endpoint[6:]
            try:
                os.unlink(socket_path)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    log.debug('could not unlink IPC endpoint {!r}: {}'.format(socket_path, e))
            else:
                log.debug('unlinked IPC endpoint {!r}'.format(socket_path))
    
    def startup(self):
        raise NotImplementedError
    
    def shutdown(self):
        raise NotImplementedError

    def __enter__(self):
        self.startup()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_traceback):
        self.shutdown()
        self._close_signal_sockets()
        try:
            self.context.close()
        except AttributeError:
            pass
        return False
    
    def _make_signal_socket(self, endpoint, socket_type=zmq.PULL):
        socket = self.context.socket(socket_type)
        socket.bind(endpoint)
        return socket
    
    def _signal_thread(self, endpoint, message='', socket_type=zmq.PUSH):
        try:
            ctlsockets = self._tls.ctlsockets
        except AttributeError:
            ctlsockets = self._tls.ctlsockets = dict()
        
        try:
            socket = ctlsockets[endpoint]
        except KeyError:
            socket = ctlsockets[endpoint] = self.context.socket(socket_type)
            socket.connect(endpoint)
                    
        socket.send(message)

    def _close_signal_sockets(self):
        try:
            ctlsockets = self._tls.__dict__.pop('ctlsockets')
        except KeyError:
            return
        
        while ctlsockets:
            _endpoint, socket = ctlsockets.popitem()
            socket.close()
            
    def __del__(self):
        self._close_signal_sockets()
        
