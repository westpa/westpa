"""
A work manager which uses ZeroMQ messaging over TCP or Unix sockets to 
distribute tasks and collect results.

The server master process streams out tasks via a PUSH socket to clients,
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

import sys, os, logging, socket, multiprocessing, threading, time, traceback, signal, random, tempfile, atexit, uuid
import cPickle as pickle
import argparse
from collections import deque
from Queue import Queue
from Queue import Empty 
import zmq
from zmq import ZMQError
try:
    from zmq.core.message import Frame
except ImportError:
    from zmq.core.message import Message as Frame

from work_managers import WorkManager, WMFuture

log = logging.getLogger(__name__)

default_ann_port     = 23811 # announcements: PUB master, SUB clients
default_task_port    = 23812 # task distribution: PUSH master, PULL clients
default_results_port = 23813 # results reception: PULL master, PUSH clients

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
        
        
        self.context = zmq.Context.instance()
                        
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
                                  .format(result.instance_id, self.instance_id))
                    
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
    def __init__(self, upstream_task_endpoint, upstream_result_endpoint, upstream_announce_endpoint,
                 nprocs = None):
        super(ZMQClient,self).__init__()
        
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_result_endpoint = upstream_result_endpoint
        self.upstream_announce_endpoint = upstream_announce_endpoint
        
        self.context = None # this really shouldn't be instantiated until after forks, just to be safe
        
        self.nprocs = nprocs or multiprocessing.cpu_count()
        
        self.workers = [] # list of worker Process objects
        self.worker_task_endpoint = None
        self.worker_result_endpoint = None
        self.active_tasks = {} # mapping of PID to active tasks
        self.worker_last_contact = {} # mapping of PID to workers
        
        self.shutdown_timeout = 1
        
        self.node_id = uuid.uuid4()
        
        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._taskfwd_ctl_endpoint = 'inproc://_taskfwd_ctl_{:x}'.format(id(self))
        self._rslfwd_ctl_endpoint  = 'inproc://_rslfwd_ctl_{:x}'.format(id(self))
        self._monitor_ctl_endpoint = 'inproc://_monitor_ctl_{:x}'.format(id(self))
        
        self._monitor_thread = None
        self._taskfwd_thread = None
        self._rslfwd_thread = None
        
        self._shutdown_signaled = False
        
    def startup(self, spawn_workers=True):
        if spawn_workers:
            self.spawn_workers()
        else:
            # primarily for testing
            self.worker_task_endpoint = self.make_ipc_endpoint()
            self.worker_result_endpoint = self.make_ipc_endpoint()            
        
        self.context = zmq.Context()
        #ctlsocket = self.context.socket(zmq.PULL)
        #ctlsocket.bind(self._startup_ctl_endpoint)
        ctlsocket = self._make_signal_socket(self._startup_ctl_endpoint)
        
        self._taskfwd_thread = threading.Thread(target=self._taskfwd_loop)
        self._taskfwd_thread.daemon = True
        self._taskfwd_thread.start()
        
        self._rslfwd_thread = threading.Thread(target=self._rslfwd_loop)
        #self._rslfwd_thread.daemon = True
        self._rslfwd_thread.start()
        
        self._monitor_thread = threading.Thread(target=self._monitor_loop)
        #self._monitor_thread.daemon = True
        self._monitor_thread.start()
        
        # Wait on all three threads starting up before returning to caller
        ctlsocket.recv()
        ctlsocket.recv()
        ctlsocket.recv()
        
        ctlsocket.close()
    
    def shutdown(self):
        if not self._shutdown_signaled:
            self._shutdown_signaled = True
            self.shutdown_workers()
            for endpoint in (self._monitor_ctl_endpoint, self._rslfwd_ctl_endpoint):
                self._signal_thread(endpoint, 'shutdown')
            
            # you may think it's a good idea to close the context here.
            # IT'S NOT -- aborts left and right
            
    
    def spawn_workers(self):
        self.worker_task_endpoint = self.make_ipc_endpoint()
        self.worker_result_endpoint = self.make_ipc_endpoint()
        for _n in xrange(self.nprocs):
            proc = ZMQWMProcess(self.worker_task_endpoint, self.worker_result_endpoint)
            self.workers.append(proc)
            proc.start()
    
    def shutdown_workers(self):
        
        # Right here is where we would implement child process pruning
        
        for proc in self.workers:
            proc.terminate()
            
        while self.workers:
            proc = self.workers.pop()
            proc.join(self.shutdown_timeout)
            if proc.is_alive():
                log.warning('sending SIGKILL to PID {}'.format(proc.pid))
                os.kill(proc.pid, signal.SIGKILL)
                proc.join()
                
            assert not proc.is_alive()
        
        
    def _monitor_loop(self):
        ctlsocket = self._make_signal_socket(self._monitor_ctl_endpoint)
        
        upstream_announce_socket = self.context.socket(zmq.SUB)
        upstream_announce_socket.setsockopt(zmq.SUBSCRIBE,'')
        upstream_announce_socket.connect(self.upstream_announce_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint)
        
        last_server_ping = 0
        
        poller = zmq.Poller()
        poller.register(upstream_announce_socket, zmq.POLLIN)
        poller.register(ctlsocket, zmq.POLLIN)
        
        try:
            while True:
                poll_results = dict(poller.poll())
                now = time.time()
                
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    if 'shutdown' in recvall(ctlsocket):
                        # this occurs during the shutdown sequence for this object, to catch
                        return
                                
                if poll_results.get(upstream_announce_socket) == zmq.POLLIN:
                    announcements = recvall(upstream_announce_socket)
                    if 'shutdown' in announcements:
                        self.shutdown()
                    elif 'ping' in announcements:
                        last_server_ping = now
        finally:
            poller.unregister(upstream_announce_socket)
            upstream_announce_socket.close(linger=0)
            
    def _taskfwd_loop(self):
        upstream_task_socket = self.context.socket(zmq.PULL)
        upstream_task_socket.connect(self.upstream_task_endpoint)
        
        worker_task_socket = self.context.socket(zmq.REP)
        worker_task_socket.bind(self.worker_task_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint)
        
        try:
            while True:
                req_frames = worker_task_socket.recv_multipart(copy=False)
                now = time.time()
                (node_id, client_id, message) = pickle.loads(req_frames[0].buffer.tobytes())
                if node_id is not None and node_id != self.node_id:
                    log.error('received request for node {}, but this is node {}'
                              .format(node_id, self.node_id))
                self.worker_last_contact[client_id] = now
                
                if message == 'task':
                    task_frames = upstream_task_socket.recv_multipart(copy=False)
                    # here is where we note what client gets what task
                    worker_task_socket.send_pyobj((self.node_id, client_id, message), flags=zmq.SNDMORE)
                    worker_task_socket.send_multipart(task_frames)
                    del task_frames
                else:
                    log.error('received invalid request from worker {!r}: {!r}'.format(client_id, message))
                    
                del req_frames                
                
        finally:
            # we never actually get here, but might as well be complete in case the code ever changes
            # to allow us to get here 
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
                now = time.time()
                
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    if 'shutdown' in recvall(ctlsocket):
                        return
                
                if poll_results.get(worker_result_socket) == zmq.POLLIN:                    
                    frames = worker_result_socket.recv_multipart(copy=False)
                    (node_id, client_id, message) = pickle.loads(frames[0].buffer.tobytes())
                    
                    
                    if node_id != self.node_id:
                        log.error('received result destined for node {} but this is node {}; ignoring'
                                  .format(node_id, self.node_id))
                        
                    self.worker_last_contact[client_id] = now                        
                        
                    if message == 'result':
                        result_meta = Result.from_zmq_frames(frames[1:], include_payload=False)
                        # here is where we note that a task has successfully finished
                        upstream_result_socket.send_multipart(frames[1:])
                        
                    else:
                        log.error('received invalid response from worker {!r}: {!r}'
                                  .format(client_id, message))
                        
                    
        finally:
            poller.unregister(worker_result_socket)
            poller.unregister(ctlsocket)
            worker_result_socket.close(linger=0)
            ctlsocket.close(linger=0)
            
        
class ZMQWorkManager(ZMQWMServer,WorkManager):
    pass


atexit.register(ZMQBase.remove_ipc_endpoints)
      
