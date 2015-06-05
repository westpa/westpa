'''
Created on May 29, 2015

@author: mzwier
'''

from __future__ import division, print_function; __metaclass__ = type

from work_managers import WMFuture

# Every ten seconds the master requests a status report from workers.
# This also notifies workers that the master is still alive
DEFAULT_STATUS_POLL = 10 

# If we haven't heard from the master or a worker (as appropriate) in these
# amounts of time, we assume a crash and shut down.
MASTER_CRASH_TIMEOUT = DEFAULT_STATUS_POLL * 6
WORKER_CRASH_TIMEOUT = DEFAULT_STATUS_POLL * 3

import logging
log = logging.getLogger(__name__)

#import gevent
import sys, uuid, socket, os,tempfile, errno, time, threading, contextlib

import signal
signames = {val:name for name, val in reversed(sorted(signal.__dict__.items()))
            if name.startswith('SIG') and not name.startswith('SIG_')}

import zmq
import numpy

def randport(address='127.0.0.1'):
    '''Select a random unused TCP port number on the given address.''' 
    s = socket.socket()
    s.bind((address,0))
    try:
        port = s.getsockname()[1]
    finally:
        s.close()
    return port

class ZMQWMError(RuntimeError):
    '''Base class for errors related to the ZeroMQ work manager itself'''
    pass

class ZMQWorkerMissing(ZMQWMError):
    '''Exception representing that a worker processing a task died or disappeared'''
    pass

class ZMQWMEnvironmentError(ZMQWMError):
    '''Class representing an error in the environment in which the ZeroMQ work manager is running.
    This includes such things as master/worker ID mismatches.'''
    
class ZMQWMTimeout(ZMQWMEnvironmentError):
    '''A timeout of a sort that indicatess that a master or worker has failed or never started.'''

class Message:
    ACK = 'ok'
    IDENTIFY = 'identify'          # Two-way identification (a reply must be an IDENTIFY message)
    TASK_REQUEST = 'task_request'
    RESULT_RETURN = 'result'
    SHUTDOWN = 'shutdown'
    MASTER_BEACON = 'master_alive'
    TASKS_AVAILABLE = 'tasks_available'
    
    VALID_MESSAGES = {ACK, IDENTIFY, 
                      TASKS_AVAILABLE, TASK_REQUEST, RESULT_RETURN,
                      SHUTDOWN, MASTER_BEACON}
    
    def __init__(self, message=None, payload=None, master_id=None, worker_id=None):
        
        if isinstance(message,Message):
            self.message   = message.message
            self.payload   = message.payload
            self.master_id = message.master_id
            self.worker_id = message.worker_id
        else:
            self.master_id = master_id
            self.worker_id = worker_id
            self.message = message
            self.payload = payload
        
    def __repr__(self):
        return ('<{!s} master_id={master_id!s} worker_id={worker_id!s} message={message!r} payload={payload!r}>'
                .format(self.__class__.__name__, **self.__dict__))   

# No need for this, since we can just put the futures in a dictionary indexed by task_id

# class TaskUpstream:
#     def __init__(self, fn, args, kwargs, future=None):
#         if future is None:
#             self.future = WMFuture()
#         else:
#             self.future = future
#             
#         self.task_id = self.future.task_id
# 
#         # Task data                    
#         self.fn = fn
#         self.args = args
#         self.kwargs = kwargs
#         
#     def __repr__(self):
#         return '<{} {task_id!s} {fn!r} {:d} args {:d} kwargs>'\
#                .format(len(self.args), len(self.kwargs), **self.__dict__)
#                
#     def __hash__(self):
#         return hash(self.task_id)
    
class Task:
    def __init__(self, fn, args, kwargs, task_id = None):
        self.task_id = task_id or uuid.uuid4()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
                
    def __repr__(self):
        return '<{} {task_id!s} {fn!r} {:d} args {:d} kwargs>'\
               .format(len(self.args), len(self.kwargs), **self.__dict__)
               
    def __hash__(self):
        return hash(self.task_id)
    
class Result:
    def __init__(self, task_id, result=None, exception=None):
        self.task_id = task_id
        self.result = result
        self.exception = exception
        
    def __repr__(self):
        return '<{} {task_id!s}>'\
               .format(**self.__dict__)
               
    def __hash__(self):
        return hash(self.task_id)
   

class PassiveTimer:
    __slots__ = {'started', 'duration'}
    def __init__(self, duration):
        self.started = time.time()
        self.duration = duration
        
    @property
    def expired(self, at=None):
        at = at or time.time()
        return (at - self.started) > self.duration
    
    @property
    def expires_in(self):
        at = time.time()
        return self.started + self.duration - at
    
    def reset(self, at=None):
        self.started = at or time.time()
        
    start = reset
    
class PassiveMultiTimer:
    def __init__(self):
        self._durations = numpy.empty((0,), float)
        self._started = numpy.empty((0,), float)
        self._identifiers = {}
        
    def add_timer(self, identifier, duration):
        
        if identifier in self._identifiers:
            raise KeyError('timer {!r} already present'.format(identifier))
        
        new_idx = len(self._durations)
        self._durations = self._durations.resize((new_idx+1,))
        self._started = self._started.resize((new_idx+1,))
        self._durations[new_idx] = duration
        self._started[new_idx] = time.time()
        
        self._identifiers[identifier] = new_idx
        
    def reset(self, identifier=None, at=None):
        at = at or time.time()
        if identifier is None:
            # reset all timers
            self._started.fill(at)
        else:
            self._started[identifier] = at
            
    
    
        

class ZMQCore:
    
    # The overall communication topology (socket layout, etc)
    # Cannot be updated without updating configuration files, command-line parameters,
    # etc. (Changes break user scripts.)
    PROTOCOL_MAJOR = 3
    
    # The set of messages and replies in use.
    # Cannot be updated without changing existing communications logic. (Changes break
    # the ZMQ WM library.)
    PROTOCOL_MINOR = 0  
    
    # Minor updates and additions to the protocol.
    # Changes do not break the ZMQ WM library, but only add new
    # functionality/code paths without changing existing code paths.    
    PROTOCOL_UPDATE = 0 
    
    PROTOCOL_VERSION = (PROTOCOL_MAJOR, PROTOCOL_MINOR, PROTOCOL_UPDATE)
    
    _ipc_endpoints_to_delete = []
    
    @classmethod    
    def make_ipc_endpoint(cls):
        (fd, socket_path) = tempfile.mkstemp()
        os.close(fd)
        endpoint = 'ipc://{}'.format(socket_path)
        cls._ipc_endpoints_to_delete.append(endpoint)
        return endpoint
    
    @classmethod
    def remove_ipc_endpoints(cls):
        while cls._ipc_endpoints_to_delete:
            endpoint = cls._ipc_endpoints_to_delete.pop()
            assert endpoint.startswith('ipc://')
            socket_path = endpoint[6:]
            try:
                os.unlink(socket_path)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    log.debug('could not unlink IPC endpoint {!r}: {}'.format(socket_path, e))
            else:
                log.debug('unlinked IPC endpoint {!r}'.format(socket_path))
    
    
    def __init__(self):
        
        # Unique identifier of this ZMQ node
        self.node_id = uuid.uuid4()
        
        # Identifier of the master node
        self.master_id = None
        
        # Beacons 
        # Master drops workers that don't check in during the beacon
        # period Workers exit if they haven't heard from master in the beacon
        # period
        self.worker_beacon_period = 5
        self.master_beacon_period = 5
        # These should allow for some fuzz, and should ratchet up as more and
        # more workers become available (maybe order 1 s for 100 workers?) This
        # should also account appropriately for startup delay on difficult
        
        
        # A friendlier description for logging
        self.node_description = '{!s} on {!s} at PID {:d}'.format(self.__class__.__name__,
                                                                  socket.gethostname(),
                                                                  os.getpid())
        
        self.validation_fail_action = 'exit' # other options are 'raise' and 'warn'
        
        self.log = logging.getLogger(__name__ + '.' + self.__class__.__name__ + '.' + str(self.node_id))
        
        # ZeroMQ context
        self.context = None
        
        # External communication endpoints
        self.rr_endpoint = None
        self.ann_endpoint = None
        
        self.inproc_endpoint = 'inproc://{!s}'.format(self.node_id)
        
        # Sockets
        self.rr_socket = None
        self.ann_socket = None
        
        # This is the main-thread end of this
        self._inproc_socket = None
        

    def __repr__(self):
        return '<{!s} {!s}>'.format(self.__class__.__name__, self.node_id)
    
    def get_identification(self):
        return {'node_id': self.node_id,
                'master_id': self.master_id,
                'class': self.__class__.__name__,
                'description': self.node_description,
                'hostname': socket.gethostname(),
                'pid': os.getpid()}
                        
    def validate_message(self, message):
        '''Validate incoming message. Raises an exception if the message is improperly formatted (TypeError)
        or does not correspond to the appropriate master (ZMQWMEnvironmentError).'''
        try:
            super_validator = super(ZMQCore,self).validate_message
        except AttributeError:
            pass
        else:
            super_validator(message)
            
        if not isinstance(message, Message):    
                raise TypeError('message is not an instance of core.Message')
        
    @contextlib.contextmanager
    def message_validation(self, msg):
        '''A context manager for message validation. The instance variable ``validation_fail_action``
        controls the behavior of this context manager:
          * 'raise': re-raise the exception that indicated failed validation. Useful for development.
          * 'exit' (default): report the error and exit the program.
          * 'warn': report the error and continue.'''
        try:
            yield
        except Exception as e:
            if self.validation_fail_action == 'raise':
                self.log.exception('message validation failed for {!r}'.format(msg))
                raise
            elif self.validation_fail_action == 'exit':
                self.log.error('message validation falied: {!s}'.format(e))
                sys.exit(1)
            elif self.validation_fail_action == 'warn':
                self.log.warn('message validation falied: {!s}'.format(e))
        
    def recv_message(self, socket, flags=0, validate=True):
        '''Receive a message object from the given socket.'''
        message = socket.recv_pyobj(flags)
        if validate:
            with self.message_validation(message):
                self.validate_message(message)
        return message
    
    def send_message(self, socket, message, payload=None, flags=0):
        '''Send a message object. Subclasses may override this to
        decorate the message with appropriate IDs, then delegate upward to actually send
        the message. ``message`` may either be a pre-constructed ``Message`` object or 
        a message identifier, in which (latter) case ``payload`` will become the message payload.
        ``payload`` is ignored if ``message`` is a ``Message`` object.'''
        
        message = Message(message, payload)
        socket.send_pyobj(message,flags)
        
    def send_reply(self, socket, original_message, reply=Message.ACK, payload=None,flags=0):
        '''Send a reply to ``original_message`` on ``socket``. The reply message
        is a Message object or a message identifier. The reply master_id and worker_id are
        set from ``original_message``, unless master_id is not set, in which case it is
        set from self.master_id.''' 
        reply = Message(reply, payload)
        reply.master_id = original_message.master_id or self.master_id
        assert original_message.worker_id is not None # should have been caught by validation prior to this
        reply.worker_id = original_message.worker_id
        self.send_message(socket, reply)
        
    def send_ack(self, socket, original_message):
        '''Send an acknowledgement message, which is mostly just to respect REQ/REP
        recv/send patterns.'''
        self.send_message(socket, Message(Message.ACK,
                                          master_id=original_message.master_id,
                                          worker_id=original_message.worker_id))

    def shutdown_handler(self, signal=None, frame=None):
        if signal is None:
            self.log.info('shutting down')
        else:
            self.log.info('shutting down on signal {!s}'.format(signames.get(signal,signal)))
        try:
            inproc_socket = self.context.socket(zmq.PUB)
            inproc_socket.connect(self.inproc_endpoint)
            self.send_message(inproc_socket, Message.SHUTDOWN, payload=signal)
        finally:
            sys.exit()
    
    def install_signal_handlers(self, signals = None):
        if not signals:
            signals = {signal.SIGINT, signal.SIGQUIT, signal.SIGTERM}
        
        for sig in signals:
            signal.signal(sig, self.shutdown_handler)

    def startup(self):
        self.context = zmq.Context()
        self.comm_thread = threading.Thread(target=self.comm_loop)
        self.comm_thread.start()
        
        self.install_signal_handlers()
      
    def shutdown(self):
        self.shutdown_handler()
        
    def join(self):
        while True:
            self.comm_thread.join(0.1)
            if not self.comm_thread.is_alive():
                break
