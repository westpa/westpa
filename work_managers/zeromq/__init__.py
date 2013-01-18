"""
A work manager which uses ZeroMQ messaging over TCP or Unix sockets to 
distribute tasks and collect results.

The server process receives task requests and results on REQ/REP sockets.
This is so that keepalive messages can be implemented if necessary. A PUB
socket is used to send critical messages -- currently, "shutdown" and "ping"
(master is alive) messages.

The client is more complex.  A client process is responsible
for starting a number of worker processes, then forwarding tasks and
results beween them and the server.  The client is also
responsible for detecting hung workers and killing them, reporting
a failure to the server. Further, each client listens for periodic
pings from the server, and will shut down if a ping is not received
in a specific time frame (indicating a crashed master).
"""

import socket
import cPickle as pickle

import zmq
from zmq import ZMQError
try:
    from zmq.core.message import Frame
except ImportError:
    from zmq.core.message import Message as Frame

#############################################
# Default values for configuration parameters
# Times are in seconds

# Time limit for tasks
# Zero means no limit
DEFAULT_TASK_TIMEOUT = 0

# When a client requests a task and none is available, wait this long to
# see if one becomes available before returning the "no task available"
# message. This is, effectively, the keepalive time for the task distribution
# connection.
DEFAULT_TASKQUEUE_WAIT = 60

# Maxmimum number of pending tasks to enqueue on server
# Zero means no limit
DEFAULT_MAX_TASKQUEUE_SIZE = 0

# How often to expect a "server is up"
DEFAULT_SERVER_HEARTBEAT_INTERVAL = 600 

# How often client wakes up to check for hung jobs or missing server
# 0 or None means to take the minimum of the server heartbeat interval
# and the task timeout.
DEFAULT_HANGCHECK_INTERVAL = 0

# How long to wait between SIGTERM and SIGKILL when terminating a worker
DEFAULT_SHUTDOWN_TIMEOUT = 5 

############
# Exceptions

class ZMQWMException(Exception):
    pass

class WorkerTerminated(ZMQWMException):
    '''Exception indicating that the worker responsible for a task was 
    terminated prior to returning results.'''
    pass


###################
# Support functions 

def randport(address='127.0.0.1'):
    '''Select a random unused port number on the given address.''' 
    s = socket.socket()
    s.bind((address,0))
    try:
        port = s.getsockname()[1]
    finally:
        s.close()
    return port

def recvall(socket):
    '''Receive and return all queued messages on the given ZeroMQ socket.'''
    messages = []
    while True:
        try:
            messages.append(socket.recv(flags=zmq.NOBLOCK))
        except ZMQError as err:
            if err.errno == zmq.EAGAIN:
                return messages
            else:
                raise

def pickle_to_frame(obj):
    '''Pickle the given object and wrap the result in a ZeroMQ
    Frame object.'''
    return Frame(pickle.dumps(obj, pickle.HIGHEST_PROTOCOL))

def unpickle_frame(frame):
    '''Unpickle the object stored in the given ZeroMQ frame in a zero-copy
    manner.'''
    if isinstance(frame, Frame):
        return pickle.loads(frame.buffer.tobytes())
    else:
        return pickle.loads(frame)
    
from . import core
from . import server, client
from . import work_manager

from server import ZMQServer
from client import ZMQClient
from work_manager import ZMQWorkManager

