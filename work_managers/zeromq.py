"""
A work manager which uses ZeroMQ messaging over TCP or Unix sockets to 
distribute tasks and collect results.

The server process streams out tasks via a PUSH socket to clients,
then receives results through a PULL socket. A PUB socket is used
to send critical messages -- currently, "shutdown" and "ping" (master is
alive) messages.

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
from collections import deque
import zmq
from zmq import ZMQError
try:
    from zmq.core.message import Frame
except ImportError:
    from zmq.core.message import Message as Frame

import work_managers
from work_managers import WorkManager, WMFuture


log = logging.getLogger(__name__)

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
    def __init__(self, server_id, task_id, fn, args, kwargs):
        self.server_id = server_id
        self.task_id = task_id
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        
    def __repr__(self):
        return '<Task {self.task_id} from server {self.server_id}: {self.fn!r}(*{self.args!r}, **{self.kwargs!r})>'\
               .format(self=self)
                
    def to_zmq_frames(self):
        header_data = {'server_id': self.server_id,
                       'task_id': self.task_id}
        payload_data = {'fn': self.fn,
                        'args': self.args,
                        'kwargs': self.kwargs}
        return [Frame(pickle.dumps(header_data, pickle.HIGHEST_PROTOCOL)),
                Frame(pickle.dumps(payload_data, pickle.HIGHEST_PROTOCOL))]
    
    @classmethod
    def from_zmq_frames(cls, frames, include_payload = True):
        assert len(frames) == 2
        header_data = pickle.loads(frames[0].buffer.tobytes())
        task = cls(header_data['server_id'], header_data['task_id'], None, None, None)
        if include_payload:
            payload_data = pickle.loads(frames[1].buffer.tobytes())
            task.fn = payload_data['fn']
            task.args = payload_data['args']
            task.kwargs = payload_data['kwargs']
            
        return task
        
class Result:
    def __init__(self, server_id, task_id, value=None, exception=None, traceback=None):
        self.server_id = server_id
        self.task_id = task_id
        self.value = value
        self.exception = exception
        self.traceback = traceback

    def __repr__(self):
        if self.exception is not None:
            return '<Result {self.task_id} for server {self.server_id}: {self.value!r})>'\
                   .format(self=self)
        else:
            return '<Result (exception) {self.task_id} for server {self.server_id}: {self.exception!r})>'\
                   .format(self=self)

        
    @property
    def is_exception(self):
        return (self.exception is not None)
    
    @property
    def is_result(self):
        return (self.exception is None)
            
    def to_zmq_frames(self):
        header_data = {'server_id': self.server_id,
                       'task_id': self.task_id,
                       'exception': self.exception,
                       'traceback': self.traceback}
        payload_data = {'value': self.value}
        return [Frame(pickle.dumps(header_data, pickle.HIGHEST_PROTOCOL)),
                Frame(pickle.dumps(payload_data, pickle.HIGHEST_PROTOCOL))]
    
    @classmethod
    def from_zmq_frames(cls, frames, include_payload = True):
        assert len(frames) == 2
        header_data = pickle.loads(frames[0].buffer.tobytes())
        result = cls(header_data['server_id'], header_data['task_id'], 
                     exception=header_data['exception'], traceback=header_data['traceback'])
        if include_payload:
            payload_data = pickle.loads(frames[1].buffer.tobytes())
            result.value = payload_data['value']
            
        return result 

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
            
class ZMQBase:
    _ipc_endpoints = []

    def __init__(self):
        # ZeroMQ context
        self.context = None
        
        # number of seconds between announcements of where to connect to the master
        self.server_heartbeat_interval = 10
        
        # This hostname
        self.hostname = socket.gethostname()
        self.host_id = '{:s}-{:d}'.format(self.hostname, os.getpid())
        self.instance_id = uuid.uuid4()

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
        return False
    
    def _make_signal_socket(self, endpoint, socket_type=zmq.PULL):
        socket = self.context.socket(socket_type)
        socket.bind(endpoint)
        return socket
    
    def _signal_thread(self, endpoint, message='', socket_type=zmq.PUSH):
        socket = self.context.socket(socket_type)
        socket.connect(endpoint)
        socket.send(message)
        socket.close()
        del socket
                
class ZMQWMServer(ZMQBase):
    
    # tasks are tuples of (task_id, function, args, keyword args)
    # results are tuples of (task_id, 'result' or 'exception', value)
    
    def __init__(self, master_task_endpoint, master_result_endpoint, master_announce_endpoint):
        super(ZMQWMServer, self).__init__()
        
        
        self.context = zmq.Context()
                        
        # where we send out work
        self.master_task_endpoint = master_task_endpoint
        
        # Where we receive results
        self.master_result_endpoint = master_result_endpoint
 
        # Where we send out announcements
        self.master_announce_endpoint = master_announce_endpoint

                        
        # tasks awaiting dispatch
        self.task_queue = deque()
        
        # futures corresponding to tasks
        self.pending_futures = dict()
        
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
        master_task_socket = self.context.socket(zmq.PUSH)
        master_task_socket.setsockopt(zmq.HWM,1)
        master_task_socket.bind(self.master_task_endpoint)

        # Create a control socket to wake up the loop        
        #ctlsocket = self.context.socket(zmq.PULL)
        #ctlsocket.bind(self._dispatch_thread_ctl_endpoint)
        ctlsocket = self._make_signal_socket(self._dispatch_thread_ctl_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint)

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        try:
            while True:
                poll_results = dict(poller.poll(100))
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if 'shutdown' in messages:
                        return
                    
                # Dispatch as many tasks as possible before checking for shutdown and waiting another .1 s
                while True:
                    try:
                        task = self.task_queue.popleft()
                    except IndexError:
                        break
                    else:
                        # this will block if no clients are around
                        master_task_socket.send_multipart(task.to_zmq_frames(), copy=False)
        finally:
            poller.unregister(ctlsocket)
            master_task_socket.close(linger=0)
            ctlsocket.close()
            
        log.debug('exiting _dispatch_loop()')
        
    def _receive_loop(self):
        # Bind the result receptor socket
        master_result_socket = self.context.socket(zmq.PULL)
        master_result_socket.bind(self.master_result_endpoint)

        # Create a control socket to wake up the loop        
        #ctlsocket = self.context.socket(zmq.PULL)
        #ctlsocket.bind(self._receive_thread_ctl_endpoint)
        ctlsocket = self._make_signal_socket(self._receive_thread_ctl_endpoint)
        
        #self._signal_startup_ctl()
        self._signal_thread(self._startup_ctl_endpoint)


        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(master_result_socket, zmq.POLLIN)
        
        try:
            while True:
                poll_results = dict(poller.poll())
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if 'shutdown' in messages:
                        return
                
                # results are tuples of (instance_id, task_id, {'result', 'exception'}, value)
                if poll_results.get(master_result_socket) == zmq.POLLIN:
                    frames = master_result_socket.recv_multipart(copy=False)
                    result = Result.from_zmq_frames(frames)
                    del frames
                    
                    if result.server_id != self.instance_id:
                        log.error('received result for instance {!s} but this is instance {!s}; ignoring. Zombie client?'
                                  .format(result.server_id, self.instance_id))
                    
                    try:
                        ft = self.pending_futures.pop(result.task_id)
                    except KeyError:
                        log.error('received result for nonexistent task {!s}; zombie client?'.format(result.task_id))
                    else:
                        if result.is_exception:
                            ft._set_exception(result.exception, result.traceback)
                        else:
                            ft._set_result(result.value)
        finally:
            poller.unregister(ctlsocket)
            poller.unregister(master_result_socket)
            master_result_socket.close(linger=0)
            ctlsocket.close()
            
        log.debug('exiting _receive_loop()')
                
    def _announce_loop(self):
        # Bind the result receptor socket
        master_announce_socket = self.context.socket(zmq.PUB)
        master_announce_socket.bind(self.master_announce_endpoint)

        # Create a control socket to wake up the loop        
        #ctlsocket = self.context.socket(zmq.PULL)
        #ctlsocket.bind(self._announce_endpoint)
        ctlsocket = self._make_signal_socket(self._announce_endpoint)
        
        #self._signal_startup_ctl()
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
                    if 'shutdown' in messages:
                        master_announce_socket.send('shutdown')
                        return
                    else:
                        for message in messages:
                            master_announce_socket.send(message)
                else:
                    # timeout
                    last_announce = time.time()
                    master_announce_socket.send('ping')
                   
                now = time.time() 
                if now - last_announce < self.server_heartbeat_interval:
                    remaining_interval = now - last_announce
                else:
                    remaining_interval = self.server_heartbeat_interval
                
        finally:
            poller.unregister(ctlsocket)
            master_announce_socket.close(linger=0)
            ctlsocket.close()
            
        log.debug('exiting _announce_loop()')
        
    def _make_append_task(self, fn, args, kwargs):
        ft = WMFuture()
        task_id = ft.task_id
        task = Task(self.instance_id, task_id, fn, args, kwargs)
        self.pending_futures[task_id] = ft
        self.task_queue.append(task)
        return ft
    
    def submit(self, fn, *args, **kwargs):
        ft = self._make_append_task(fn, args, kwargs)
        # wake up the dispatch loop
        self._signal_thread(self._dispatch_thread_ctl_endpoint)
        return ft
    
    def submit_many(self, tasks):
        futures = [self._make_append_task(fn, args, kwargs) for (fn,args,kwargs) in tasks]
        self._signal_thread(self._dispatch_thread_ctl_endpoint)
        return futures
        
    def shutdown(self):
        if not self._shutdown_signaled:
            self._shutdown_signaled = True
            
            for endpoint in (self._dispatch_thread_ctl_endpoint,self._receive_thread_ctl_endpoint,self._announce_endpoint):
                self._signal_thread(endpoint, 'shutdown')            

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
        
        result_socket = self.context.socket(zmq.PUSH)
        result_socket.setsockopt(zmq.HWM,1) # block, rather than queue, when waiting to dispatch results
        result_socket.connect(self.upstream_result_endpoint)
                        
        associated_server_id = None
        associated_node_id = None
        this_client_id = os.getpid()
                
        try:
            while True:
                
                # Get a task with request/reply to identify ourselves with which task we get, for
                # future bookkeeping
                req_tuple = (associated_node_id, this_client_id, 'task')
                task_socket.send_pyobj(req_tuple)
                
                # do a zero-copy receive and unpickling to minimize RAM use
                reply_frames = task_socket.recv_multipart(copy=False)
                (returned_node_id, returned_client_id, message) = pickle.loads(reply_frames[0].buffer.tobytes())

                
                # Make sure we're talking to the correct client process
                if associated_node_id is None:
                    associated_node_id = returned_node_id
                elif returned_node_id != associated_node_id:
                    raise ValueError('received reply from node {} when expecting reply from {}'
                                     .format(returned_node_id, associated_node_id))
                
                # Make sure that the client process is talking to us
                if returned_client_id != this_client_id:
                    raise ValueError('received reply destined for {} (this is client process {})'
                                     .format(returned_client_id, this_client_id))
                
                if message == 'task':
                    task = Task.from_zmq_frames(reply_frames[1:])
                    
                    if associated_server_id is None:
                        associated_server_id = task.server_id
                    elif task.server_id != associated_server_id:
                        raise ValueError('received task from server {} when expecting task from server {}'
                                         .format(task.server_id, associated_server_id))
                    
                                            
                    try:
                        result_value = task.fn(*task.args, **task.kwargs)
                    except Exception as e:
                        result_object = Result(task.server_id, task.task_id, exception=e, traceback=traceback.format_exc())
                    else:
                        result_object = Result(task.server_id, task.task_id, value=result_value)
    
                        
                    # submit a result
                    
                    result_socket.send_pyobj((associated_node_id, this_client_id, 'result'), flags=zmq.SNDMORE)
                    result_socket.send_multipart(result_object.to_zmq_frames())
                    
                else:
                    log.error('unknown message {!r} received'.format(message))
                    
                del reply_frames
                    
        finally:
            task_socket.close(linger=0)
            result_socket.close(linger=0)

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
                              help='''Use the given ZeroMQ endpoint for task distribution. (Use {argname} over
                                      explicit endpoints, if possible.)'''.format(argname=wmenv.arg_flag('zmq_comm_mode')))
        wm_group.add_argument(wmenv.arg_flag('zmq_result_endpoint'), metavar='RESULT_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for result collection. (Use {argname} over
                                      explicit endpoints, if possible.)'''.format(argname=wmenv.arg_flag('zmq_comm_mode')))
        wm_group.add_argument(wmenv.arg_flag('zmq_announce_endpoint'), metavar='ANNOUNCE_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for task distribution. (Use {argname} over
                                      explicit endpoints, if possible.)'''.format(argname=wmenv.arg_flag('zmq_comm_mode')))
        wm_group.add_argument(wmenv.arg_flag('zmq_task_timeout'), metavar='TIMEOUT', type=int,
                              help='''Kill worker processes that take longer than TIMEOUT seconds''')
                                             
    
    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 
        
        n_workers = wmenv.get_val('n_workers', multiprocessing.cpu_count(), int)
        hangcheck = wmenv.get_val('zmq_task_timeout', 60, int)


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
        
        return cls(task_endpoint, result_endpoint, announce_endpoint, n_workers, hangcheck)
                     
    
    def __init__(self, upstream_task_endpoint, upstream_result_endpoint, upstream_announce_endpoint,
                 n_workers = None, hangcheck=60):
        super(ZMQClient,self).__init__()
        
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
        self.shutdown_timeout = 1
        
        # How long we wait for a response from a worker process before
        # declaring it hung and terminating it -- None means do not check.
        self.worker_task_timeout = None
        
        # How often (in s) we check for a hung worker
        self.hangcheck = hangcheck
                
        self.node_id = uuid.uuid4()
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
            self.running = True   
            if spawn_workers:
                with self.worker_lock:
                    for _n in xrange(self.n_workers):
                        self._spawn_worker()
                
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
                self._signal_thread(endpoint, 'shutdown')
                    
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
    
            # these are all atomic operations
            assert proc.pid in self.workers
            del self.workers[proc.pid]
            active_task = self.worker_active_tasks.pop(proc.pid,None)
            self.worker_last_dispatch.pop(proc.pid,None)
            
            if active_task:
                upstream_result_socket = self.context.socket(zmq.PUSH)
                upstream_result_socket.connect(self.upstream_result_endpoint)
                
                try:
                    result = Result(active_task.server_id, active_task.task_id,
                                    exception=WorkerTerminated('worker performing this task was terminated'),
                                    traceback='')

                    upstream_result_socket.send_multipart(result.to_zmq_frames())
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
                poll_results = dict(poller.poll(min(self.hangcheck,self.server_heartbeat_interval)*1000))
                now = time.time()
                
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if 'shutdown' in messages:
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
                    if 'shutdown' in announcements:
                        self._shutdown()
                    elif 'ping' in announcements:
                        last_server_ping = now
                
                if last_server_ping is not None and (now-last_server_ping) > self.server_heartbeat_interval:
                    log.error('no communication from server; shutting down')
                    self._shutdown()
                    
                if self.worker_task_timeout:
                    # this is not atomic, so borderline cases undergoing updates may trigger
                    # false positive hangs...but only borderline cases, which is sufficient
                    # for us
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
            
    def _taskfwd_loop(self):
        ctlsocket = self._make_signal_socket(self._taskfwd_ctl_endpoint)
        
        upstream_task_socket = self.context.socket(zmq.PULL)
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
        
        try:
            while True:
                worker_poll_status = dict(worker_poller.poll())
                
                if worker_poll_status.get(ctlsocket) == zmq.POLLIN:
                    if 'shutdown' in recvall(ctlsocket):
                        return
                    
                if worker_poll_status.get(worker_task_socket) == zmq.POLLIN:
                    req_frames = worker_task_socket.recv_multipart(copy=False)
                    (node_id, client_id, message) = pickle.loads(req_frames[0].buffer.tobytes())
                    if node_id is not None and node_id != self.node_id:
                        log.error('received request for node {}, but this is node {}'
                                  .format(node_id, self.node_id))
                    
                    if message == 'task':
                        upstream_poll_status = dict(upstream_poller.poll())
                        
                        if upstream_poll_status.get(ctlsocket) == zmq.POLLIN:
                            if 'shutdown' in recvall(ctlsocket):
                                return
                            
                        assert upstream_poll_status.get(upstream_task_socket) == zmq.POLLIN
                        
                        task_frames = upstream_task_socket.recv_multipart(copy=False)
                        task_meta = Task.from_zmq_frames(task_frames, include_payload=False)
                                                
                        worker_task_socket.send_pyobj((self.node_id, client_id, message), flags=zmq.SNDMORE)
                        worker_task_socket.send_multipart(task_frames)
                        
                        self.worker_active_tasks[client_id] = task_meta
                        self.worker_last_dispatch[client_id] = time.time()
                        
                        del task_frames
                    else:
                        log.error('received invalid request from worker {!r}: {!r}'.format(client_id, message))
                        
                    del req_frames                
                
        finally:
            upstream_poller.unregister(ctlsocket)
            upstream_poller.unregister(upstream_task_socket)
            worker_poller.unregister(ctlsocket)
            worker_poller.unregister(worker_task_socket)
            worker_task_socket.close(linger=0)
            upstream_task_socket.close(linger=0)
            
    def _rslfwd_loop(self):
        ctlsocket = self._make_signal_socket(self._rslfwd_ctl_endpoint)
        upstream_result_socket = self.context.socket(zmq.PUSH)
        upstream_result_socket.setsockopt(zmq.HWM, 1)
        upstream_result_socket.connect(self.upstream_result_endpoint)
        worker_result_socket = self.context.socket(zmq.PULL)
        worker_result_socket.bind(self.worker_result_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint)
        
        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(worker_result_socket, zmq.POLLIN)
        
        try:
            while True:
                poll_results = dict(poller.poll())
                
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    if 'shutdown' in recvall(ctlsocket):
                        return
                
                if poll_results.get(worker_result_socket) == zmq.POLLIN:                    
                    frames = worker_result_socket.recv_multipart(copy=False)
                    (node_id, client_id, message) = pickle.loads(frames[0].buffer.tobytes())
                    
                    if node_id != self.node_id:
                        log.error('received result destined for node {} but this is node {}; ignoring'
                                  .format(node_id, self.node_id))                            
                    if message == 'result':
                        result_meta = Result.from_zmq_frames(frames[1:], include_payload=False)                        
                        upstream_result_socket.send_multipart(frames[1:])
                    else:
                        log.error('received invalid response from worker {!r}: {!r}'
                                  .format(client_id, message))
                        
                    
        finally:
            poller.unregister(worker_result_socket)
            poller.unregister(ctlsocket)
            worker_result_socket.close(linger=0)
            ctlsocket.close(linger=0)

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
        wm_group.add_argument(wmenv.arg_flag('zmq_comm_mode'), metavar='COMM_MODE', choices=('tcp', 'ipc'),
                              help='''Use TCP/IP ({argname}=tcp) or Unix ({argname}=ipc) sockets for communication.
                                    IPC sockets are more efficient for single-node communication, but do not allow 
                                    communication between nodes. Default is TCP.'''.format(argname=wmenv.arg_flag('zmq_comm_mode')))
        wm_group.add_argument(wmenv.arg_flag('zmq_task_endpoint'), metavar='TASK_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for task distribution. (Use {argname} over
                                      explicit endpoints, if possible.)'''.format(argname=wmenv.arg_flag('zmq_comm_mode')))
        wm_group.add_argument(wmenv.arg_flag('zmq_result_endpoint'), metavar='RESULT_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for result collection. (Use {argname} over
                                      explicit endpoints, if possible.)'''.format(argname=wmenv.arg_flag('zmq_comm_mode')))
        wm_group.add_argument(wmenv.arg_flag('zmq_announce_endpoint'), metavar='ANNOUNCE_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for task distribution. (Use {argname} over
                                      explicit endpoints, if possible.)'''.format(argname=wmenv.arg_flag('zmq_comm_mode')))
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
        server_info_filename = wmenv.get_val('zmq_server_info', 'zmq_server_info_{}.json'.format(uuid.uuid4().hex))
        
        # if individual endpoints are named, we use these
        tests = [not bool(wmenv.get_val('zmq_task_endpoint')),
                 not bool(wmenv.get_val('zmq_result_endpoint')),
                 not bool(wmenv.get_val('zmq_announce_endpoint'))]
        if all(tests):
            # No endpoints specified; see if we have been instructed to choose TCP or IPC
            comm_mode = wmenv.get_val('zmq_comm_mode')
            if not comm_mode: comm_mode = 'tcp'
            comm_mode = comm_mode.lower()
            if comm_mode not in ('tcp', 'ipc'):
                raise ValueError('invalid ZMQ communications mode: {!r}'.format(comm_mode))
            elif comm_mode == 'tcp':
                # Choose random ports
                task_endpoint = cls.canonicalize_endpoint('tcp://*')
                result_endpoint = cls.canonicalize_endpoint('tcp://*')
                announce_endpoint = cls.canonicalize_endpoint('tcp://*')
            else: # ipc
                task_endpoint = cls.make_ipc_endpoint()
                result_endpoint = cls.make_ipc_endpoint()
                announce_endpoint = cls.make_ipc_endpoint()
        elif any(tests):
            raise ValueError('either none or all three endpoints must be specified')
        else:
            task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_task_endpoint'))
            result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_result_endpoint'))
            announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_announce_endpoint'))
            
        return cls(n_workers, task_endpoint, result_endpoint, announce_endpoint, server_info_filename = server_info_filename,
                   hangcheck=hangcheck)
        
    def remove_server_info_file(self):
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
                 write_server_info = True, server_info_filename=None, hangcheck=60):
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

        ZMQWMServer.__init__(self, master_task_endpoint, master_result_endpoint, master_announce_endpoint)            
        
        if n_workers > 0:
            # this node is both a master and a client; start workers
            self.internal_client = ZMQClient(master_task_endpoint, master_result_endpoint, master_announce_endpoint,
                                             self.n_workers, hangcheck)
        
        if write_server_info:
            self.server_info_filename = server_info_filename or 'zmq_server_info_{}.json'.format(uuid.uuid4().hex)
            self.write_server_info(self.server_info_filename)
            atexit.register(self.remove_server_info_file)
                
                    
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
      
