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


from __future__ import division, print_function; __metaclass__ = type

import os, logging, socket, multiprocessing, threading, time, traceback, signal, tempfile, atexit, uuid, json, re, random
import cPickle as pickle
import Queue as queue
import zmq
from zmq import ZMQError
try:
    from zmq.core.message import Frame
except ImportError:
    from zmq.core.message import Message as Frame

import work_managers
from work_managers import WorkManager, WMFuture

log = logging.getLogger(__name__)

# General message format: tuples of the form (tag, server_id, client_id, payload),
# where tag is a string, server_id is a UUID or None, client_id is a UUID or 
# None, and payload is any associated data, or None. This applies to each frame.

# Control messages
MSG_ACK = 'ok'
MSG_SHUTDOWN = 'shutdown'
MSG_PING = 'ping' # ping inquiry
MSG_PONG = 'pong' # ping reply
MSG_TASK_REQUEST = 'task'
MSG_TASK_AVAILABLE = 'task_avail'
MSG_TASK_UNAVAILABLE = 'task_unavail'
MSG_RESULT_SUBMISSION = 'result'

RESULT_TYPE_RETVAL = 'retval'
RESULT_TYPE_EXCEPTION = 'exception'


def randport():
    s = socket.socket()
    s.bind(('127.0.0.1',0))
    port = s.getsockname()[1]
    s.close()
    return port

class ZMQWMException(Exception):
    pass

class WorkerTerminated(ZMQWMException):
    '''Exception indicating that the worker responsible for a task was 
    terminated prior to returning results.'''
    pass

class Task:
    def __init__(self, fn, args, kwargs, future=None):
        if future is None:
            self.future = WMFuture()
        else:
            self.future = future
            
        self.task_id = self.future.task_id

        # Task data                    
        self.fn = fn
        self.args = args
        self.kwargs = kwargs                

def recvall(socket):
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
    
ShutdownSentinel = object()

            
class ZMQBase:
    _ipc_endpoints = []

    def __init__(self, server_heartbeat_interval = 10):
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
            
                        
class ZMQWMServer(ZMQBase):
    
    def __init__(self, master_task_endpoint, master_result_endpoint, master_announce_endpoint,
                 server_heartbeat_interval, max_taskqueue_size=None, taskqueue_wait=10):
        super(ZMQWMServer, self).__init__(server_heartbeat_interval)
        
        self.context = zmq.Context()
                        
        # where we send out work
        self.master_task_endpoint = master_task_endpoint
        
        # Where we receive results
        self.master_result_endpoint = master_result_endpoint
 
        # Where we send out announcements
        self.master_announce_endpoint = master_announce_endpoint
        
        # tasks awaiting dispatch
        self.task_queue = queue.Queue(max_taskqueue_size or 0)
        self.task_queue_wait = taskqueue_wait
        
        # tasks awaiting return, keyed by ID
        self.pending_tasks = dict()
        
        self._shutdown_signaled = False

        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._dispatch_thread_ctl_endpoint = 'inproc://_dispatch_thread_ctl_{:x}'.format(id(self))        
        self._receive_thread_ctl_endpoint = 'inproc://_receive_thread_ctl_{:x}'.format(id(self))
        self._announce_endpoint = 'inproc://_announce_{:x}'.format(id(self))

    def write_server_info(self, filename):
        pass
        
    def startup(self):
        # start up server threads, blocking until their sockets are ready
        
        # create an inproc socket to sequence the startup of worker threads
        # each thread needs to write an empty message to this endpoint so
        # that startup() doesn't exit until all required sockets are open
        # and listening
        
        #ctlsocket = self.context.socket(zmq.PULL)
        #ctlsocket.bind(self._startup_ctl_endpoint)
        ctlsocket = self._make_signal_socket(self._startup_ctl_endpoint)
        
        #proper use here is to start a thread, then recv, and in the thread func
        #use _signal_startup_ctl() once all its sockets are bound
        self._dispatch_thread = threading.Thread(target=self._dispatch_loop)
        self._dispatch_thread.start()        
        
        self._receive_thread = threading.Thread(target=self._receive_loop)
        self._receive_thread.start()
        
        self._announce_thread = threading.Thread(target=self._announce_loop)
        self._announce_thread.start()
        
        ctlsocket.recv() # dispatch
        ctlsocket.recv() # receive
        ctlsocket.recv() # announce
        
        ctlsocket.close()
        
                        
    def _dispatch_loop(self):
        # a 1-1 mapping between items added to the task queue and ZMQ sends
        
        # Bind the task distributor socket
        master_task_socket = self.context.socket(zmq.REP)
        master_task_socket.bind(self.master_task_endpoint)

        # Create a control socket to wake up the loop for shutdown 
        ctlsocket = self._make_signal_socket(self._dispatch_thread_ctl_endpoint)
        
        # Signal that this thread has started
        self._signal_thread(self._startup_ctl_endpoint)

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(master_task_socket, zmq.POLLIN)
        
        debug_logging = log.isEnabledFor(logging.DEBUG)
        info_logging = log.isEnabledFor(logging.INFO)
        this_server_id = self.instance_id
        
        try:
            while True:                    
                poll_results=dict(poller.poll())
                
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if MSG_SHUTDOWN in messages:
                        return
                
                if poll_results.get(master_task_socket) == zmq.POLLIN:
                    message = master_task_socket.recv_pyobj()
                    if debug_logging:
                        log.debug('server: received message {!r}'.format(message))
                    (tag, server_id, client_id, _payload) = message[:4]
                    if server_id and server_id != this_server_id:
                        log.error('received message {!r} destined for another server ({!s}); rogue/zombie client?'
                                  .format(tag, server_id))
                    
                    if tag == MSG_TASK_REQUEST:
                        # dispatch task, or 'no task available' message
                        try:
                            # Block up to one second
                            # TODO: make this a configurable parameter
                            task = self.task_queue.get(block=True, timeout=self.task_queue_wait)
                        except queue.Empty:
                            log.debug('server: no task available')
                            master_task_socket.send_pyobj((MSG_TASK_UNAVAILABLE, this_server_id, client_id, None))
                        else:
                            if task is ShutdownSentinel:
                                return
                            
                            # Send task as a two-part message: a (small) metadata header
                            # and a (possibly large) payload
                            log.debug('server: dispatching task {!s} ({!r})'.format(task.task_id, task.fn))
                            master_task_socket.send_pyobj((MSG_TASK_AVAILABLE, this_server_id, client_id, task.task_id),
                                                          flags=zmq.SNDMORE)
                            master_task_socket.send(pickle_to_frame((task.fn, task.args, task.kwargs)), copy=False)
                    else:
                        log.error('unknown/unsupported message received on task socket: {!r}'.format(tag))
        finally:
            poller.unregister(ctlsocket)
            poller.unregister(master_task_socket)
            master_task_socket.close(linger=0)
            ctlsocket.close()
            log.debug('server: exiting dispatch loop')
        
    def _receive_loop(self):
        # Bind the result receptor socket
        master_result_socket = self.context.socket(zmq.REP)
        master_result_socket.bind(self.master_result_endpoint)

        # Create a control socket to wake up the loop        
        ctlsocket = self._make_signal_socket(self._receive_thread_ctl_endpoint)

        # Signal that this thread has started        
        self._signal_thread(self._startup_ctl_endpoint)

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(master_result_socket, zmq.POLLIN)

        debug_logging = log.isEnabledFor(logging.DEBUG)
        info_logging = log.isEnabledFor(logging.INFO)
        this_server_id = self.instance_id
        
        try:
            while True:
                poll_results = dict(poller.poll())
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if MSG_SHUTDOWN in messages:
                        return
                
                # results are tuples of (instance_id, task_id, 'result' or 'exception', value)
                if poll_results.get(master_result_socket) == zmq.POLLIN:
                    frames = master_result_socket.recv_multipart(copy=False)
                    message = unpickle_frame(frames[0])
                    if debug_logging:
                        log.debug('server: received message {!r}'.format(message))
                        
                    (tag, server_id, client_id, payload) = message[:4]
                    if server_id and server_id != this_server_id:
                        log.error('received message {!r} destined for another server ({!s}); rogue/zombie client?'
                                  .format(tag, server_id))
                        
                    master_result_socket.send_pyobj((MSG_ACK, this_server_id, client_id, None))

                    if tag == MSG_RESULT_SUBMISSION:
                        task_id = payload
                        try:
                            task = self.pending_tasks[task_id]
                        except KeyError:
                            log.error('received result for unknown task {!s}'.format(task_id))
                        else:
                            assert task.task_id == task_id
                            (result_type, result_payload) = unpickle_frame(frames[1])
                            if result_type == RESULT_TYPE_EXCEPTION:
                                (exception, traceback) = result_payload
                                if debug_logging:
                                    log.debug('server: received exception for task {!s} ({!r}): {!s}'
                                              .format(task_id, task.fn, exception))
                                task.future._set_exception(exception, traceback)
                                del exception, traceback
                            elif result_type == RESULT_TYPE_RETVAL:
                                retval = result_payload
                                if debug_logging:
                                    log.debug('server: received result for task {!s} ({!r})'.format(task_id, retval))
                                task.future._set_result(retval)
                            else:
                                log.error('unknown result type received for task {!s} ({!r})'.format(task_id, task.fn))
                            del result_payload
                    else:
                        log.error('unknown/unsupported message received on result socket: {!r}'.format(tag))
        finally:
            poller.unregister(ctlsocket)
            poller.unregister(master_result_socket)
            master_result_socket.close(linger=0)
            ctlsocket.close()
            log.debug('server: exiting receive loop')
                
    def _announce_loop(self):
        # Bind the result receptor socket
        master_announce_socket = self.context.socket(zmq.PUB)
        master_announce_socket.bind(self.master_announce_endpoint)

        # Create a control socket to wake up the loop        
        ctlsocket = self._make_signal_socket(self._announce_endpoint)
        
        # Signal that this thread has started
        self._signal_thread(self._startup_ctl_endpoint)

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        
        last_announce = 0
        remaining_interval = self.server_heartbeat_interval
        
        try:
            while True:
                poll_results = dict(poller.poll(remaining_interval*1000))
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if MSG_SHUTDOWN in messages:
                        log.debug('server: shutdown message received; forwarding to clients')
                        master_announce_socket.send(MSG_SHUTDOWN)
                        return
                    else:
                        for message in messages:
                            master_announce_socket.send(message)
                            
                now = time.time()
                if now - last_announce >= self.server_heartbeat_interval:
                    last_announce = now
                    master_announce_socket.send(MSG_PING)
                    remaining_interval = self.server_heartbeat_interval
                else:
                    remaining_interval = self.server_heartbeat_interval - (now - last_announce)                            
        finally:
            poller.unregister(ctlsocket)
            master_announce_socket.close(linger=0)
            ctlsocket.close()
            log.debug('server: exiting announce loop')
            
    def submit(self, fn, args=None, kwargs=None):
        return self.submit_many([(fn,args if args is not None else (),kwargs if kwargs is not None else {})])[0]
    
    def submit_many(self, tasks):
        futures = []
        pending_tasks = self.pending_tasks
        task_queue = self.task_queue
        
        for (fn,args,kwargs) in tasks:
            task = Task(fn, args, kwargs)
            pending_tasks[task.task_id] = task
            futures.append(task.future)
            task_queue.put(task)            
        return futures
        
    def shutdown(self):
        if not self._shutdown_signaled:
            self._shutdown_signaled = True
            
            for endpoint in (self._dispatch_thread_ctl_endpoint,self._receive_thread_ctl_endpoint,self._announce_endpoint):
                self._signal_thread(endpoint, MSG_SHUTDOWN)
                
            # Put a sentinel on the task queue to wake up waits on it
            self.task_queue.put(ShutdownSentinel)
                        

class ZMQWMProcess(ZMQBase,multiprocessing.Process):
    '''A worker process, meant to be run via multiprocessing.Process()'''
    
    def __init__(self, upstream_task_endpoint, upstream_result_endpoint):
        ZMQBase.__init__(self)
        multiprocessing.Process.__init__(self)   
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_result_endpoint = upstream_result_endpoint
                
    def run(self):
        '''Run a recieve work/do work/dispatch result loop. This is *designed* to hang
        in the event of a hung task function, as the parent process is responsible
        for managing the worker process pool, forcefully if necessary.'''
        
        # block SIGINT, so that we can only shut down if told so by an announcement, or
        # if we are sent SIGQUIT. This avoids some hangs when running on only one node
        # and the user presses CTRL-C
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        
        self.context = zmq.Context()
        
        task_socket = self.context.socket(zmq.REQ)
        task_socket.connect(self.upstream_task_endpoint)
        
        result_socket = self.context.socket(zmq.REQ)
        result_socket.setsockopt(zmq.HWM,1) # block, rather than queue, when waiting to dispatch results
        result_socket.connect(self.upstream_result_endpoint)
                        
        associated_server_id = None
        associated_node_id = None
        this_client_id = self.instance_id
        this_client_pid = os.getpid()

        debug_logging = log.isEnabledFor(logging.DEBUG)
        info_logging = log.isEnabledFor(logging.INFO)
                
        try:
            while True:
                inner_req_tuple = (MSG_TASK_REQUEST, associated_server_id, this_client_id, None)
                outer_req_tuple = (MSG_TASK_REQUEST, associated_node_id, this_client_id, this_client_pid)
                task_socket.send_pyobj(outer_req_tuple, flags=zmq.SNDMORE)
                task_socket.send_pyobj(inner_req_tuple)
                
                # do a zero-copy receive and unpickling to minimize RAM use
                task_frames = task_socket.recv_multipart(copy=False)
                
                # task_frames[0]: node wrapper
                # task_frames[1]: server header
                # task_frames[2]: task payload (if any)
                
                node_header = unpickle_frame(task_frames[0])
                if debug_logging:
                    log.debug('process: node header: {!r}'.format(node_header))

                (_, node_id, _, client_pid) = node_header
                
                # Make sure we're talking to the correct client process
                if associated_node_id is None:
                    associated_node_id = node_id
                elif node_id != associated_node_id:
                    raise ValueError('received reply from node {} when expecting reply from {}'.format(node_id, associated_node_id))
                
                # Make sure that the client process is talking to whom it thinks it should
                if client_pid != this_client_pid:
                    raise ValueError('received reply destined for PID {} (this is client process {})'
                                     .format(client_pid, this_client_pid))
                    
                # Unpack header
                message = unpickle_frame(task_frames[1])
                if debug_logging:
                    log.debug('process: message: {!r}'.format(message))
                (tag, server_id, client_id, task_id) = message
                
                # Check server ID
                if associated_server_id is None:
                    associated_server_id = server_id
                elif server_id != associated_server_id:
                    raise ValueError('received response from server {!s} when expecting task from server {!s}'
                                     .format(server_id, associated_server_id))
                
                # Check that the server is talking to whom it thinks it should
                if client_id != this_client_id:
                    raise ValueError('received response for client {!s}, but this is client {!s}'
                                     .format(client_id, this_client_id))
                
                if tag == MSG_TASK_AVAILABLE:
                    fn, args, kwargs = unpickle_frame(task_frames[2])
                elif tag == MSG_TASK_UNAVAILABLE:
                    log.debug('process: no task available')
                    continue
                else:
                    raise ValueError('unknown message {!r} received'.format(tag))
                
                del node_header, message, task_frames
                
                # We only get here if we are error-free and actually have a task to run
                log.debug('process: running task {!s} ({!r})'.format(task_id, fn))
                try:
                    result_type = RESULT_TYPE_RETVAL
                    result_payload = fn(*args,**kwargs)
                except Exception as e:
                    result_type = RESULT_TYPE_EXCEPTION
                    result_payload = (e, traceback.format_exc())
                    log.debug('process: task {!s} ({!r}) failed with exception: {!s}'.format(task_id, fn, e))
                else:
                    log.debug('process: task {!s} ({!r}) completed successfully'.format(task_id, fn))
                    
                result_socket.send_pyobj((MSG_RESULT_SUBMISSION, associated_node_id, this_client_id, this_client_pid, task_id),
                                         flags=zmq.SNDMORE)
                result_socket.send_pyobj((MSG_RESULT_SUBMISSION, associated_server_id, this_client_id, task_id), flags=zmq.SNDMORE)
                result_socket.send_pyobj((result_type, result_payload))
                
                # Receive acknowledgement
                result_socket.recv()
                    
        finally:
            task_socket.close(linger=0)
            result_socket.close(linger=0)
            log.debug('process: exiting run loop')

class ZMQClient(ZMQBase):
    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 

        wm_group = parser.add_argument_group('options for ZeroMQ ("zmq") client')
        wm_group.add_argument(wmenv.arg_flag('zmq_server_info'), metavar='SERVER_INFO_FILE',
                              help='Store server information (if master) or obtain server information (if client) '
                                   'from SERVER_INFO_FILE. This is helpful if running server and clients on multiple '
                                   'machines which share a filesystem, as explicit hostnames/ports are not required')
        wm_group.add_argument(wmenv.arg_flag('zmq_task_endpoint'), metavar='TASK_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for task distribution.''')
        wm_group.add_argument(wmenv.arg_flag('zmq_result_endpoint'), metavar='RESULT_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for result collection.''')
        wm_group.add_argument(wmenv.arg_flag('zmq_announce_endpoint'), metavar='ANNOUNCE_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for task distribution.''')
        wm_group.add_argument(wmenv.arg_flag('zmq_heartbeat_interval'), metavar='INTERVAL',
                              help='''If a client has not
                                      heard from the server in approximately INTERVAL seconds, the client will
                                      assume the server has crashed and shut down. This may need to be increased
                                      from the default on heavily loaded systems. (Default: 60 seconds.)''')
        wm_group.add_argument(wmenv.arg_flag('zmq_task_timeout'), metavar='TIMEOUT', type=int,
                              help='''Kill worker processes that take longer than TIMEOUT seconds''')
                                             
    
    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 
        
        n_workers = wmenv.get_val('n_workers', multiprocessing.cpu_count(), int)
        hangcheck = wmenv.get_val('zmq_task_timeout', 60, int)
        heartbeat = wmenv.get_val('zmq_heartbeat_interval', 60, int)


        tests = [not bool(wmenv.get_val('zmq_task_endpoint')),
                 not bool(wmenv.get_val('zmq_result_endpoint')),
                 not bool(wmenv.get_val('zmq_announce_endpoint'))]
        if all(tests):
            # No endpoints specified; use server info file
            server_info_filename = wmenv.get_val('zmq_server_info')
            if server_info_filename is None:
                raise EnvironmentError('neither endpoints nor server info file specified')
            else:
                try:
                    server_info = json.load(open(server_info_filename,'rt'))
                    task_endpoint = server_info['task_endpoint']
                    result_endpoint = server_info['result_endpoint']
                    announce_endpoint = server_info['announce_endpoint']    
                except Exception as e:
                    raise EnvironmentError('cannot load server info file {!r}: {}'.format(server_info_filename,e))                
        elif any(tests):
            raise ValueError('either none or all three endpoints must be specified')
        else:
            task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_task_endpoint'),allow_wildcard_host=False)
            result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_result_endpoint'),allow_wildcard_host=False)
            announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_announce_endpoint'),allow_wildcard_host=False)
        
        return cls(task_endpoint, result_endpoint, announce_endpoint, n_workers, hangcheck, heartbeat)
                     
    
    def __init__(self, upstream_task_endpoint, upstream_result_endpoint, upstream_announce_endpoint,
                 n_workers = None, task_timeout=60, server_heartbeat_interval = 10):
        super(ZMQClient,self).__init__(server_heartbeat_interval)
        
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_result_endpoint = upstream_result_endpoint
        self.upstream_announce_endpoint = upstream_announce_endpoint
        
        self.context = None # this really shouldn't be instantiated until after forks, just to be safe
        
        self.n_workers = n_workers or multiprocessing.cpu_count()
        
        self.worker_task_endpoint = self.make_ipc_endpoint()
        self.worker_result_endpoint = self.make_ipc_endpoint()
                
        self.worker_lock = threading.RLock()
        self.workers = {} # mapping of PID to Process objects
        self.worker_active_tasks = {} # mapping of PID to Task objects (without payloads)
        self.worker_last_dispatch = {} # mapping of PID to time when last dispatch occurred
        
        # How long we wait for worker processes to shutdown on SIGTERM
        # before moving on with SIGKILL
        self.shutdown_timeout = 10
        
        # How long we wait for a response from a worker process before
        # declaring it hung and terminating it -- None means do not check.
        self.worker_task_timeout = task_timeout
        
        # How often (in s) we check for a hung worker
        self.hangcheck = 30
                
        self.associated_server_id = None
        
        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._taskfwd_ctl_endpoint = 'inproc://_taskfwd_ctl_{:x}'.format(id(self))
        self._rslfwd_ctl_endpoint  = 'inproc://_rslfwd_ctl_{:x}'.format(id(self))
        self._monitor_ctl_endpoint = 'inproc://_monitor_ctl_{:x}'.format(id(self))
        
        self._monitor_thread = None
        self._taskfwd_thread = None
        self._rslfwd_thread = None
        
        self._shutdown_signaled = False
        
        self.running = False
        
    def startup(self, spawn_workers=True):     
        if not self.running:
            from work_managers import environment
            self.running = True   
            if spawn_workers:
                with self.worker_lock:
                    pi_name = '{}_PROCESS_INDEX'.format(environment.WMEnvironment.env_prefix)
                    for n in xrange(self.n_workers):
                        os.environ[pi_name] = str(n)
                        self._spawn_worker()
                    try:
                        del os.environ[pi_name]
                    except KeyError:
                        pass
                
            self.context = zmq.Context()
            ctlsocket = self._make_signal_socket(self._startup_ctl_endpoint)
            
            try:
                self._taskfwd_thread = threading.Thread(target=self._taskfwd_loop)
                self._taskfwd_thread.start()
                
                self._rslfwd_thread = threading.Thread(target=self._rslfwd_loop)
                self._rslfwd_thread.start()
                
                self._monitor_thread = threading.Thread(target=self._monitor_loop)
                self._monitor_thread.start()
                
                # Wait on all three threads starting up before continuing
                ctlsocket.recv()
                ctlsocket.recv()
                ctlsocket.recv()
    
            finally:
                ctlsocket.close()
            
    def _shutdown(self):
        if not self._shutdown_signaled:
            self._shutdown_signaled = True
            for endpoint in (self._monitor_ctl_endpoint, self._rslfwd_ctl_endpoint, self._taskfwd_ctl_endpoint):
                self._signal_thread(endpoint, MSG_SHUTDOWN)
                    
    def shutdown(self):
        if self.running:
            self._shutdown()
            self._wait_for_shutdown()
            self.running = False
        
    def run(self):
        self._wait_for_shutdown()
            
    def _wait_for_shutdown(self):
        self._monitor_thread.join()
        self._rslfwd_thread.join()
        self._taskfwd_thread.join()
        
    def _spawn_worker(self):
        with self.worker_lock:
            proc = ZMQWMProcess(self.worker_task_endpoint, self.worker_result_endpoint)
            proc.start()
            assert proc.pid is not None
            self.workers[proc.pid] = proc        
    
    def _shutdown_all_workers(self):
        
        # Right here is where we would implement zombie child process slaying
        # e.g. using psutil package
        
        worker_procs = self.workers.values()
        shutdown_threads = []
        
        for proc in worker_procs:
            thread = threading.Thread(target=self._shutdown_worker, args=(proc,))
            thread.start()
            shutdown_threads.append(thread)
            
        for thread in shutdown_threads:
            thread.join()
            
        assert not self.workers
                
    def _shutdown_worker(self, proc):
        with self.worker_lock:
            assert proc.is_alive()
            proc.terminate()
            proc.join(self.shutdown_timeout)
            if proc.is_alive():
                log.warning('sending SIGKILL to PID {}'.format(proc.pid))
                os.kill(proc.pid, signal.SIGKILL)
                proc.join()
            assert not proc.is_alive()
    
            assert proc.pid in self.workers
            del self.workers[proc.pid]
            active_task_id = self.worker_active_tasks.pop(proc.pid,None)
            self.worker_last_dispatch.pop(proc.pid,None)
            
            if active_task_id is not None:
                upstream_result_socket = self.context.socket(zmq.REQ)
                upstream_result_socket.connect(self.upstream_result_endpoint)
                
                try:
                    upstream_result_socket.send_pyobj((MSG_RESULT_SUBMISSION, None, self.instance_id, active_task_id),
                                                      flags=zmq.SNDMORE)
                    upstream_result_socket.send_pyobj((RESULT_TYPE_EXCEPTION,
                                                       (WorkerTerminated('worker performing this task was terminated'), '')))
                finally:
                    upstream_result_socket.close()
            
        
    def _monitor_loop(self):
        ctlsocket = self._make_signal_socket(self._monitor_ctl_endpoint)
        
        upstream_announce_socket = self.context.socket(zmq.SUB)
        upstream_announce_socket.setsockopt(zmq.SUBSCRIBE,'')
        upstream_announce_socket.connect(self.upstream_announce_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint)
        
        last_server_ping = None
        
        poller = zmq.Poller()
        poller.register(upstream_announce_socket, zmq.POLLIN)
        poller.register(ctlsocket, zmq.POLLIN)
        
        try:
            while True:
                waittime = min(self.hangcheck,self.server_heartbeat_interval)*1000
                poll_results = dict(poller.poll(waittime))
                now = time.time()
                
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if MSG_SHUTDOWN in messages:
                        self._shutdown_all_workers()
                        return
                    # other messages are directives to manage workers
                    for message in messages:
                        fields = message.split(':')
                        if fields[0] == 'terminate':
                            self._shutdown_worker(self.workers[int(fields[1])])
                        elif fields[1] == 'spawn':
                            self._spawn_worker()
                                
                if poll_results.get(upstream_announce_socket) == zmq.POLLIN:
                    announcements = recvall(upstream_announce_socket)
                    if MSG_SHUTDOWN in announcements:
                        log.debug('client: shutdown received')
                        self._shutdown()
                    elif MSG_PING in announcements:
                        log.debug('client: ping received')
                        last_server_ping = now
                
                if last_server_ping is not None and (now-last_server_ping) > 3*self.server_heartbeat_interval:
                    log.error('no communication from server; shutting down')
                    self._shutdown()
                    
                if self.worker_task_timeout:
                    with self.worker_lock:
                        last_dispatches = self.worker_last_dispatch.copy()
                        for (pid,last_dispatch) in last_dispatches.iteritems():
                            if (now - last_dispatch) > self.worker_task_timeout:
                                log.error('worker at PID {} has hung; terminating and starting new worker'
                                          .format(pid))
                                self._signal_thread(self._monitor_ctl_endpoint, 'terminate:{}'.format(pid))
                                self._signal_thread(self._monitor_ctl_endpoint, 'spawn:')                            
                                        
        finally:
            poller.unregister(upstream_announce_socket)
            upstream_announce_socket.close(linger=0)
            log.debug('client: exiting monitor loop')
            
    def _taskfwd_loop(self):
        ctlsocket = self._make_signal_socket(self._taskfwd_ctl_endpoint)
        
        upstream_task_socket = self.context.socket(zmq.REQ)
        upstream_task_socket.connect(self.upstream_task_endpoint)
        
        worker_task_socket = self.context.socket(zmq.REP)
        worker_task_socket.bind(self.worker_task_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint)
        
        worker_poller = zmq.Poller()
        worker_poller.register(ctlsocket, zmq.POLLIN)
        worker_poller.register(worker_task_socket, zmq.POLLIN)
        
        upstream_poller = zmq.Poller()
        upstream_poller.register(ctlsocket, zmq.POLLIN)
        upstream_poller.register(upstream_task_socket, zmq.POLLIN)

        debug_logging = log.isEnabledFor(logging.DEBUG)
        info_logging = log.isEnabledFor(logging.INFO)
        
        this_node_id = self.instance_id
        
        try:
            while True:
                worker_poll_status = dict(worker_poller.poll())
                
                if worker_poll_status.get(ctlsocket) == zmq.POLLIN:
                    if MSG_SHUTDOWN in recvall(ctlsocket):
                        return
                    
                if worker_poll_status.get(worker_task_socket) == zmq.POLLIN:
                    downstream_req_frames = worker_task_socket.recv_multipart(copy=False)
                    
                    #inner_req_tuple = (MSG_TASK_REQUEST, associated_server_id, this_client_id, None)
                    #outer_req_tuple = (MSG_TASK_REQUEST, associated_node_id, this_client_id, this_client_pid)
                    
                    node_header = unpickle_frame(downstream_req_frames[0])
                    if debug_logging:
                        log.debug('node header: {!r}'.format(node_header))    
                    (downstream_tag, node_id, client_id, client_pid) = node_header
                    
                    if node_id is not None and node_id != this_node_id:
                        raise ValueError('received request for node {}, but this is node {}'
                                         .format(node_id, this_node_id))
                                            
                    if downstream_tag == MSG_TASK_REQUEST:
                        while True:
                            upstream_task_socket.send_multipart(downstream_req_frames[1:])
                            upstream_poll_status = dict(upstream_poller.poll())
                            
                            if upstream_poll_status.get(ctlsocket) == zmq.POLLIN:
                                if MSG_SHUTDOWN in recvall(ctlsocket):
                                    return
                            elif upstream_poll_status.get(upstream_task_socket) == zmq.POLLIN:
                                upstream_response_frames = upstream_task_socket.recv_multipart(copy=False)
                                upstream_message = unpickle_frame(upstream_response_frames[0])
                                if debug_logging:
                                    log.debug('received message from upstream: {!r}'.format(upstream_message))
                                (upstream_tag, _, _, payload) = upstream_message
                                
                                if upstream_tag == MSG_TASK_AVAILABLE:
                                    task_id = payload
                                    worker_task_socket.send_pyobj((upstream_tag, this_node_id, client_id, client_pid),
                                                                  flags=zmq.SNDMORE)
                                    worker_task_socket.send_multipart(upstream_response_frames)

                                    with self.worker_lock:
                                        self.worker_active_tasks[client_pid] = task_id
                                        self.worker_last_dispatch[client_pid] = time.time()
                                    break # from inner loop
                                elif upstream_message == MSG_TASK_UNAVAILABLE:
                                    log.debug('no task available')
                                    continue #inner loop
                    else:
                        raise ValueError('received invalid request from worker {!s} (PID {!s}): {!r}'
                                         .format(client_id, client_pid, downstream_tag))
        finally:
            upstream_poller.unregister(ctlsocket)
            upstream_poller.unregister(upstream_task_socket)
            worker_poller.unregister(ctlsocket)
            worker_poller.unregister(worker_task_socket)
            worker_task_socket.close(linger=0)
            upstream_task_socket.close(linger=0)
            log.debug('exiting task fowarding loop')
            
    def _rslfwd_loop(self):
        ctlsocket = self._make_signal_socket(self._rslfwd_ctl_endpoint)
        upstream_result_socket = self.context.socket(zmq.REQ)
        upstream_result_socket.connect(self.upstream_result_endpoint)
        worker_result_socket = self.context.socket(zmq.REP)
        worker_result_socket.bind(self.worker_result_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint)
        
        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(worker_result_socket, zmq.POLLIN)

        debug_logging = log.isEnabledFor(logging.DEBUG)
        info_logging = log.isEnabledFor(logging.INFO)
        
        this_node_id = self.instance_id
        
        try:
            while True:
                poll_results = dict(poller.poll())
                
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    if MSG_SHUTDOWN in recvall(ctlsocket):
                        return
                elif poll_results.get(worker_result_socket) == zmq.POLLIN:
                    downstream_frames = worker_result_socket.recv_multipart(copy=False)
                    worker_result_socket.send_pyobj((MSG_ACK,None,None,None))
                    node_header = unpickle_frame(downstream_frames[0])
                    if debug_logging:
                        log.debug('node header: {!r}'.format(node_header))    
                    (downstream_tag, node_id, client_id, client_pid, task_id) = node_header
                    
                    if node_id is not None and node_id != this_node_id:
                        raise ValueError('received request for node {}, but this is node {}'
                                         .format(node_id, this_node_id))
                    
                    if downstream_tag == MSG_RESULT_SUBMISSION:
                        # remove outstanding task ID for this client
                        with self.worker_lock:
                            try:
                                active_task_id = self.worker_active_tasks.pop(client_pid)
                            except KeyError:
                                # may have been replaced by newer task
                                log.debug('pid {} not found in active tasks {!r}'.format(client_pid,
                                                                                         self.worker_active_tasks))
                            else:
                                if active_task_id != task_id:
                                    self.worker_active_tasks[client_pid] = active_task_id
                            
                        upstream_result_socket.send_multipart(downstream_frames[1:])
                        upstream_result_socket.recv()
                    else:
                        raise ValueError('received invalid response from worker {!r}: {!r}'
                                         .format(client_pid, downstream_tag))
        finally:
            poller.unregister(worker_result_socket)
            poller.unregister(ctlsocket)
            worker_result_socket.close(linger=0)
            ctlsocket.close(linger=0)
            log.debug('exiting result forwarding loop')

    @property
    def is_master(self):
        '''True if this is the master process for task distribution. This is necessary, e.g., for
        MPI, where all processes start identically and then must branch depending on rank.'''
        return False
            

class ZMQWorkManager(ZMQWMServer,WorkManager):
    write_server_info = True
    
    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 

        wm_group = parser.add_argument_group('options for ZeroMQ ("zmq") work manager')
        wm_group.add_argument(wmenv.arg_flag('zmq_mode'), metavar='MODE', choices=('server', 'client'),
                              help='Operate as a server (MODE=server) or a client (MODE=client).')
        wm_group.add_argument(wmenv.arg_flag('zmq_server_info'), metavar='SERVER_INFO_FILE',
                              help='Store server information (if master) or obtain server information (if client) '
                                   'from SERVER_INFO_FILE. This is helpful if running server and clients on multiple '
                                   'machines which share a filesystem, as explicit hostnames/ports are not required')
        wm_group.add_argument(wmenv.arg_flag('zmq_task_endpoint'), metavar='TASK_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for task distribution.''')
        wm_group.add_argument(wmenv.arg_flag('zmq_result_endpoint'), metavar='RESULT_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for result collection.''')
        wm_group.add_argument(wmenv.arg_flag('zmq_announce_endpoint'), metavar='ANNOUNCE_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for task distribution.''')
        wm_group.add_argument(wmenv.arg_flag('zmq_heartbeat_interval'), metavar='INTERVAL',
                              help='''If a client has not
                                      heard from the server in approximately INTERVAL seconds, the client will
                                      assume the server has crashed and shut down. This may need to be increased
                                      from the default on heavily loaded systems. (Default: 10 seconds.)''')
        wm_group.add_argument(wmenv.arg_flag('zmq_task_timeout'), metavar='TIMEOUT', type=int,
                              help='''Kill worker processes that take longer than TIMEOUT seconds''')
            
            
    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 
         
        if wmenv.get_val('zmq_mode','server').lower() == 'client':
            return ZMQClient.from_environ()
        
        n_workers = wmenv.get_val('n_workers', multiprocessing.cpu_count(), int)
        hangcheck = wmenv.get_val('zmq_task_timeout', 60, int)
        heartbeat = wmenv.get_val('zmq_heartbeat_interval', 10, int)
        server_info_filename = wmenv.get_val('zmq_server_info', 'zmq_server_info_{}.json'.format(uuid.uuid4().hex))
        
        # if individual endpoints are named, we use these
        tests = [not bool(wmenv.get_val('zmq_task_endpoint')),
                 not bool(wmenv.get_val('zmq_result_endpoint')),
                 not bool(wmenv.get_val('zmq_announce_endpoint'))]
        if all(tests):
            # Choose random ports
            task_endpoint = cls.canonicalize_endpoint('tcp://*')
            result_endpoint = cls.canonicalize_endpoint('tcp://*')
            announce_endpoint = cls.canonicalize_endpoint('tcp://*')
        elif any(tests):
            raise ValueError('either none or all three endpoints must be specified')
        else:
            task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_task_endpoint'))
            result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_result_endpoint'))
            announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_announce_endpoint'))
            
        return cls(n_workers, task_endpoint, result_endpoint, announce_endpoint, server_info_filename = server_info_filename,
                   hangcheck=hangcheck, server_heartbeat_interval=heartbeat)
        
    def remove_server_info_file(self):
        if self.server_info_filename:
            filename = self.server_info_filename
            try:
                os.unlink(filename)
            except OSError as e:
                log.debug('could not remove server info file {!r}: {}'.format(filename, e))
            else:
                log.debug('removed server info file {!r}'.format(filename))
    
    def write_server_info(self, filename=None):
        filename = filename or self.server_info_filename
        hostname = socket.gethostname()
        with open(filename, 'wt') as infofile:
            json.dump({'task_endpoint': re.sub(r'\*', hostname, self.master_task_endpoint),
                       'result_endpoint': re.sub(r'\*', hostname, self.master_result_endpoint),
                       'announce_endpoint': re.sub(r'\*', hostname, self.master_announce_endpoint)},
                      infofile)
        os.chmod(filename, 0600)
    
    def __init__(self, n_workers = None, 
                 master_task_endpoint = None, master_result_endpoint = None, master_announce_endpoint = None,
                 write_server_info = True, server_info_filename=None, hangcheck=60,
                 server_heartbeat_interval = 10):
        WorkManager.__init__(self)
        
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()
        self.n_workers = n_workers
        self.internal_client = None
        self.hangcheck = hangcheck
        
        argtests = [master_task_endpoint is None, master_result_endpoint is None, master_announce_endpoint is None]
        if any(argtests) and not all(argtests):
            raise ValueError('endpoints must all be either specified or None (not mixed)')
        else:
            assign_endpoints = all(argtests)
            
        if assign_endpoints:
            master_task_endpoint = self.make_ipc_endpoint()
            master_result_endpoint = self.make_ipc_endpoint()
            master_announce_endpoint = self.make_ipc_endpoint()

        ZMQWMServer.__init__(self, master_task_endpoint, master_result_endpoint, master_announce_endpoint,
                             server_heartbeat_interval)            
        
        if n_workers > 0:
            # this node is both a master and a client; start workers
            self.internal_client = ZMQClient(master_task_endpoint, master_result_endpoint, master_announce_endpoint,
                                             self.n_workers, hangcheck)
        
        if write_server_info:
            self.server_info_filename = server_info_filename or 'zmq_server_info_{}.json'.format(uuid.uuid4().hex)
            self.write_server_info(self.server_info_filename)
            atexit.register(self.remove_server_info_file)
        else:
            self.server_info_filename = None
                
                    
    def startup(self):
        if not self.running:
            self.running = True
            ZMQWMServer.startup(self)
            if self.internal_client is not None:
                self.internal_client.startup()
            
    def shutdown(self):
        if self.running:
            if self.internal_client is not None:
                self.internal_client.shutdown()
            ZMQWMServer.shutdown(self)
            self.remove_server_info_file()
            self.running = False


atexit.register(ZMQBase.remove_ipc_endpoints)
      
