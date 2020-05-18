'''
Created on May 29, 2015

@author: mzwier
'''

from pickle import UnpicklingError

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
import sys, uuid, socket, os,tempfile, errno, time, threading, contextlib, traceback, multiprocessing, json, re
from collections import OrderedDict

import signal
signames = {val:name for name, val in reversed(sorted(signal.__dict__.items()))
            if name.startswith('SIG') and not name.startswith('SIG_')}

import zmq
import numpy

DEFAULT_LINGER = 1

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
    SHUTDOWN = 'shutdown'
    
    ACK = 'ok'
    NAK = 'no'
    IDENTIFY = 'identify'          # Two-way identification (a reply must be an IDENTIFY message)
    TASKS_AVAILABLE = 'tasks_available'
    TASK_REQUEST = 'task_request'
    
    
    MASTER_BEACON = 'master_alive'
    RECONFIGURE_TIMEOUT = 'reconfigure_timeout'    
    
    TASK = 'task'
    RESULT = 'result'
    
    idempotent_announcement_messages = {SHUTDOWN, TASKS_AVAILABLE, MASTER_BEACON}

    
    def __init__(self, message=None, payload=None, master_id=None, src_id=None):
        
        if isinstance(message,Message):
            self.message    = message.message
            self.payload    = message.payload
            self.master_id  = message.master_id
            self.src_id     = message.src_id
        else:
            self.master_id  = master_id
            self.src_id     = src_id
            self.message = message
            self.payload = payload
        
    def __repr__(self):
        return ('<{!s} master_id={master_id!s} src_id={src_id!s} message={message!r} payload={payload!r}>'
                .format(self.__class__.__name__, **self.__dict__))
        
        
    @classmethod
    def coalesce_announcements(cls, messages):
        d = OrderedDict()
        for msg in messages:
            if msg.message in cls.idempotent_announcement_messages:
                key = msg.message
            else:
                key = (msg.message, msg.payload)
            d[key] = msg
        coalesced = list(msg.values())
        log.debug('coalesced {} announcements into {}'.format(len(messages), len(coalesced)))
        return coalesced

TIMEOUT_MASTER_BEACON = 'master_beacon'
TIMEOUT_WORKER_CONTACT = 'worker_contact'
               
class Task:
    def __init__(self, fn, args, kwargs, task_id = None):
        self.task_id = task_id or uuid.uuid4()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
                
    def __repr__(self):
        try:
            return '<{} {task_id!s} {fn!r} {:d} args {:d} kwargs>'\
                   .format(self.__class__.__name__, len(self.args), len(self.kwargs), **self.__dict__)
        except TypeError:
            # no length
            return '<{} {task_id!s} {fn!r}'.format(self.__class__.__name__, **self.__dict__)
               
    def __hash__(self):
        return hash(self.task_id)
    
    def execute(self):
        '''Run this task, returning a Result object.'''
        rsl = Result(task_id = self.task_id)
        try:
            rsl.result = self.fn(*self.args, **self.kwargs)
        except BaseException as e:
            rsl.exception = e
            rsl.traceback = traceback.format_exc()
        return rsl
        
    
class Result:
    def __init__(self, task_id, result=None, exception=None, traceback=None):
        self.task_id = task_id
        self.result = result
        self.exception = exception
        self.traceback = traceback
        
    def __repr__(self):
        return '<{} {task_id!s} ({})>'\
               .format(self.__class__.__name__, 'result' if self.exception is None else 'exception', **self.__dict__)
               
    def __hash__(self):
        return hash(self.task_id)
   

class PassiveTimer:
    __slots__ = {'started', 'duration'}
    def __init__(self, duration, started=None):
        if started is None:
            started = time.time()
        self.started = started
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
        self._identifiers = numpy.empty((0,), numpy.object_)
        self._durations = numpy.empty((0,), float)
        self._started = numpy.empty((0,), float)
        self._indices = {} # indexes into durations/started, keyed by identifier
        
    def add_timer(self, identifier, duration):
        
        if identifier in self._identifiers:
            raise KeyError('timer {!r} already present'.format(identifier))
        
        new_idx = len(self._identifiers)
        self._durations.resize((new_idx+1,))
        self._started.resize((new_idx+1,))
        self._identifiers.resize((new_idx+1,))
        self._durations[new_idx] = duration
        self._started[new_idx] = time.time()
        self._identifiers[new_idx] = identifier
        self._indices[identifier] = new_idx
        
    def remove_timer(self, identifier):
        idx = self._indices.pop(identifier)
        self._durations = numpy.delete(self._durations, idx)
        self._started = numpy.delete(self._started, idx)
        self._identifiers = numpy.delete(self._identifiers, idx)
        
    def change_duration(self, identifier, duration):
        idx = self._indices[identifier]
        self._durations[idx] = duration
        
    def reset(self, identifier=None, at=None):
        at = at or time.time()
        if identifier is None:
            # reset all timers
            self._started.fill(at)
        else:
            self._started[self._indices[identifier]] = at
    
    def expired(self, identifier, at = None):
        at = at or time.time()
        idx = self._indices[identifier]
        return (at - self._started[idx]) > self._durations[idx]
    
    def next_expiration(self):
        at = time.time()
        idx = (self._started + self._durations - at).argmin()
        return self._identifiers[idx]
    
    def next_expiration_in(self):
        at = time.time()
        idx = (self._started + self._durations - at).argmin()
        next_at = self._started[idx] + self._durations[idx] - at
        return next_at if next_at > 0 else 0
    
    def which_expired(self, at=None):
        at = at or time.time()
        expired_indices = (at - self._started) > self._durations
        return self._identifiers[expired_indices]
        
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
    
    # The default transport for "internal" (inter-thread/-process) communication
    # IPC should work except on really odd systems with no local storage
    internal_transport = 'ipc'
    
    default_comm_mode = 'ipc'
    default_master_heartbeat = 20.0
    default_worker_heartbeat = 20.0
    default_timeout_factor = 5.0
    default_startup_timeout = 120.0
    default_shutdown_timeout = 5.0
        
    
    
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
                
    @classmethod
    def make_tcp_endpoint(cls, address='127.0.0.1'):
        return 'tcp://{}:{}'.format(address,randport(address))
    
    @classmethod
    def make_internal_endpoint(cls):
        assert cls.internal_transport in {'ipc', 'tcp'}
        if cls.internal_transport == 'ipc':
            return cls.make_ipc_endpoint()
        else: # cls.internal_transport == 'tcp'
            return cls.make_tcp_endpoint()
    
    def __init__(self):
        
        # Unique identifier of this ZMQ node
        self.node_id = uuid.uuid4()
        
        # Identifier of the task distribution network (work manager)
        self.network_id = None
        
        # Beacons
        # Workers expect to hear from the master at least every master_beacon_period
        # Master expects to hear from the workers at least every worker_beacon_period
        # If more than {master,worker}_beacon_period*timeout_factor elapses, the 
        # master/worker is considered missing.
         
        self.worker_beacon_period = self.default_worker_heartbeat
        self.master_beacon_period = self.default_master_heartbeat
        self.timeout_factor = self.default_timeout_factor
        
        # These should allow for some fuzz, and should ratchet up as more and
        # more workers become available (maybe order 1 s for 100 workers?) This
        # should also account appropriately for startup delay on difficult
        # systems.
        
        # Number of seconds to allow first contact between at least one worker
        # and the master.
        self.startup_timeout = self.default_startup_timeout
        
        
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
        
        self.master_id = None
        
        if os.environ.get('WWMGR_ZMQ_DEBUG_MESSAGES', 'n').upper() in {'Y', 'YES', '1', 'T', 'TRUE'}:
            self._super_debug = True
        else:
            self._super_debug = None

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
        if message.src_id is None:
            raise ZMQWMEnvironmentError('message src_id is not set')
        if self.master_id is not None and message.master_id is not None and message.master_id != self.master_id:
            raise ZMQWMEnvironmentError('incoming message associated with another master (this={!s}, incoming={!s}'.format(self.master_id, message.master_id))
                
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
    
    def recv_message(self, socket, flags=0, validate=True, timeout=None):
        '''Receive a message object from the given socket, using the given flags.
        Message validation is performed if ``validate`` is true.
        If ``timeout`` is given, then it is the number of milliseconds to wait
        prior to raising a ZMQWMTimeout exception. ``timeout`` is ignored if
        ``flags`` includes ``zmq.NOBLOCK``.'''
        
        if timeout is None or flags & zmq.NOBLOCK:
            message = socket.recv_pyobj(flags)
        else:        
            poller = zmq.Poller()
            poller.register(socket, zmq.POLLIN)
            try:
                poll_results = dict(poller.poll(timeout=timeout))
                if socket in poll_results:
                    message = socket.recv_pyobj(flags)
                else:
                    raise ZMQWMTimeout('recv timed out')
            finally:
                poller.unregister(socket)
        
        if self._super_debug:
            self.log.debug('received {!r}'.format(message))
        if validate:
            with self.message_validation(message):
                self.validate_message(message)
        return message
    
    def recv_all(self, socket, flags=0, validate=True):
        '''Receive all messages currently available from the given socket.'''
        messages = []
        while True:
            try:
                messages.append(self.recv_message(socket, flags | zmq.NOBLOCK, validate))
            except zmq.Again:
                return messages
            
    def recv_ack(self, socket, flags=0, validate=True, timeout=None):
        msg = self.recv_message(socket, flags, validate, timeout)
        if validate:
            with self.message_validation(msg):
                assert msg.message in (Message.ACK, Message.NAK)
        return msg
    
    def send_message(self, socket, message, payload=None, flags=0):
        '''Send a message object. Subclasses may override this to
        decorate the message with appropriate IDs, then delegate upward to actually send
        the message. ``message`` may either be a pre-constructed ``Message`` object or 
        a message identifier, in which (latter) case ``payload`` will become the message payload.
        ``payload`` is ignored if ``message`` is a ``Message`` object.'''
        
        message = Message(message, payload)
        if message.master_id is None:
            message.master_id = self.master_id
        message.src_id=self.node_id
        
        if self._super_debug:
            self.log.debug('sending {!r}'.format(message))
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
                                          master_id=original_message.master_id or self.master_id,
                                          src_id=self.node_id))
        
    def send_nak(self, socket, original_message):
        '''Send a negative acknowledgement message.'''
        self.send_message(socket, Message(Message.NAK,
                                          master_id=original_message.master_id or self.master_id,
                                          src_id=self.node_id))

    def send_inproc_message(self, message, payload=None, flags=0):
        inproc_socket = self.context.socket(zmq.PUB)
        inproc_socket.connect(self.inproc_endpoint)
        # annoying wait for sockets to settle
        time.sleep(0.01)
        self.send_message(inproc_socket, message, payload, flags)
        # used to be a close with linger here, but it was cutting off messages
        
    def signal_shutdown(self):
        try:
            self.send_inproc_message(Message.SHUTDOWN)
        except AttributeError:
            # this is expected if self.context has been set to None (i.e. it has already been destroyed)
            pass
        except Exception as e:
            self.log.debug('ignoring exception {!r} in signal_shutdown()'.format(e))

    def shutdown_handler(self, signal=None, frame=None):
        if signal is None:
            self.log.info('shutting down')
        else:
            self.log.info('shutting down on signal {!s}'.format(signames.get(signal,signal)))
        self.signal_shutdown()
    
    def install_signal_handlers(self, signals = None):
        if not signals:
            signals = {signal.SIGINT, signal.SIGQUIT, signal.SIGTERM}
        
        for sig in signals:
            signal.signal(sig, self.shutdown_handler)
            
    def install_sigint_handler(self):
        self.install_signal_handlers()

    def startup(self):
        self.context = zmq.Context()
        self.comm_thread = threading.Thread(target=self.comm_loop)
        self.comm_thread.start()
        
        #self.install_signal_handlers()
      
    def shutdown(self):
        self.shutdown_handler()
        
    def join(self):
        while True:
            self.comm_thread.join(0.1)
            if not self.comm_thread.is_alive():
                break


def shutdown_process(process, timeout=1.0):
    process.join(timeout)
    if process.is_alive():            
        log.debug('sending SIGINT to process {:d}'.format(process.pid))
        os.kill(process.pid, signal.SIGINT)
        process.join(timeout)
        if process.is_alive():
            log.warning('sending SIGKILL to worker process {:d}'.format(process.pid))
            os.kill(process.pid, signal.SIGKILL)
            process.join()
            
        log.debug('process {:d} terminated with code {:d}'.format(process.pid, process.exitcode))
    else:
        log.debug('worker process {:d} terminated gracefully with code {:d}'.format(process.pid, process.exitcode))        
    assert not process.is_alive()    

class IsNode:
    def __init__(self, n_local_workers=None):
        from work_managers.zeromq.worker import ZMQWorker
        
        if n_local_workers is None:
            n_local_workers = multiprocessing.cpu_count()
        
        self.downstream_rr_endpoint = None
        self.downstream_ann_endpoint = None

        if n_local_workers:
            self.local_ann_endpoint = self.make_internal_endpoint()
            self.local_rr_endpoint = self.make_internal_endpoint()
            self.local_workers = [ZMQWorker(self.local_rr_endpoint, self.local_ann_endpoint) for _n in range(n_local_workers)]
        else:
            self.local_ann_endpoint = None
            self.local_rr_endpoint = None
            self.local_workers = []
            
        self.local_worker_processes = [multiprocessing.Process(target = worker.startup, args=(n,)) 
                                       for (n, worker) in enumerate(self.local_workers)]
        
        self.host_info_files = []           

    def write_host_info(self, filename=None):
        filename = filename or 'zmq_host_info_{}.json'.format(self.node_id.hex)
        hostname = socket.gethostname()
                        
        with open(filename, 'wt') as infofile:
            info = {}
            info['rr_endpoint'] = re.sub(r'\*', hostname, self.downstream_rr_endpoint or '')
            info['ann_endpoint'] = re.sub(r'\*', hostname, self.downstream_ann_endpoint or '')
            json.dump(info,infofile)
        self.host_info_files.append(filename)

    def startup(self):
        for process in self.local_worker_processes:
            process.start()
            
    def shutdown(self):
        try:
            shutdown_timeout = self.shutdown_timeout
        except AttributeError:
            shutdown_timeout = 1.0
            
        for process in self.local_worker_processes:
            shutdown_process(process, shutdown_timeout)
            
        for host_info_file in self.host_info_files:
            try:
                os.unlink(host_info_file)
            except OSError:
                pass
            
