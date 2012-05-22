"""
A work manager which uses ZeroMQ messaging over TCP or Unix sockets to 
distribute tasks and collect results.

The server master process streams out tasks via a PUSH socket to clients,
then receives results through a PULL socket. A PUB socket is used
to send critical messages -- currently, "shutdown" and "ping" (master is
alive) messages.

The client process is more complex.  A "master client" is responsible
for starting a number of worker processes, then forwarding tasks and
results between them and the server.  The master client is also
responsible for detecting hung workers and killing them, reporting
a failure to the server. Further, each client listens for periodic
pings from the server, and will shut down if a ping is not received
in a specific time frame (indicating a crashed master).

Task messages are pickled tuples of the form::

  (instance_id, task_id, fn, args, kwargs)
  
where ``instance_id`` is an identifier indicating the server which has
dispatched the instance (for detecting duplicate servers), ``task_id`` is
an ID uniquely identifying the task, and the remaining elements are the
task itself, which is executed as ``fn(*args,**kwargs)``.

Results messages are pickled tuples of the form::

  (instance_id, task_id, result_type, payload)

where ``instance_id`` must be the instance identifier of the originating
server process (to check for crashed and restarted master processes),
``task_id`` is the task to which the result corresponts, ``result_type``
is either 'result' indicating a successful invocation of the task 
function, or 'exception' indicating that an exception occurred executing the
task function. ``payload`` is either the return value of the task function,
or an ``(exception, traceback_str)`` tuple, where ``exception`` is the 
exception raised while calling the task function, and ``traceback_str`` is
a formatted traceback.

Announcement messages are simple strings, either 'shutdown' to direct clients
to shut down, or 'ping' to announce that the server is still alive.

"""


from __future__ import division, print_function; __metaclass__ = type

import sys, os, logging, socket, multiprocessing, threading, time, traceback, signal, random, tempfile, atexit, uuid
import argparse
from collections import deque
from Queue import Queue
from Queue import Empty 
import zmq
from zmq import ZMQError

from work_managers import WorkManager, WMFuture

log = logging.getLogger(__name__)

default_ann_port     = 23811 # announcements: PUB master, SUB clients
default_task_port    = 23812 # task distribution: PUSH master, PULL clients
default_results_port = 23813 # results reception: PULL master, PUSH clients

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
                
    def get_port(self, host='127.0.0.1'):
        '''Get an unused port. This is not atomic, but may be sufficient in certain 
        cases.'''
        s = socket.socket()
        s.bind((host,0))
        port = s.getsockname()[1]
        s.close()
        return port
    
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
    
    def _signal_thread(self, endpoint, message='', socket_type=zmq.PUSH):
        socket = self.context.socket(socket_type)
        socket.connect(endpoint)
        socket.send(message)
        socket.close()
        del socket
    
class ZMQWMServer(ZMQBase):
    
    # tasks are tuples of (task_id, function, args, keyword args)
    # results are tuples of (task_id, {'result', 'exception'}, value)
    
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
        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(self._startup_ctl_endpoint)
        
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
        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(self._dispatch_thread_ctl_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint, socket_type=zmq.PUSH)

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        try:
            while True:
                poll_results = dict(poller.poll(100))
                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if 'shutdown' in messages:
                        return
                    
                # Run as many tasks as possible before checking for shutdown and waiting another .1 s
                while True:
                    try:
                        task_tuple = self.task_queue.popleft()
                    except IndexError:
                        break
                    else:
                        # this will block if no clients are around
                        master_task_socket.send_pyobj(task_tuple)
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
        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(self._receive_thread_ctl_endpoint)
        
        #self._signal_startup_ctl()
        self._signal_thread(self._startup_ctl_endpoint, socket_type=zmq.PUSH)


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
                    try:
                        (instance_id, task_id, result_type, payload) = master_result_socket.recv_pyobj()
                    except ValueError:
                        log.error('received malformed result; ignorning')
                    else:
                        if instance_id != self.instance_id:
                            log.error('received result for instance {!s} but this is instance {!s}; zombie process?'
                                      .format(instance_id, self.instance_id))
                        
                        try:
                            ft = self.pending_futures.pop(task_id)
                        except KeyError:
                            log.error('received result for nonexistent task {!s}; zombie client?'.format(task_id))
                        else: 
                            if result_type == 'result':
                                ft._set_result(payload)
                            elif result_type == 'exception':
                                ft._set_exception(*payload)
                            else:
                                log.error('received unknown result type {!r} for task {!s}; incompatible/zombie client?'
                                          .format(result_type, task_id))
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
        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(self._announce_endpoint)
        
        #self._signal_startup_ctl()
        self._signal_thread(self._startup_ctl_endpoint, socket_type=zmq.PUSH)

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
        task_tuple = (self.instance_id, task_id, fn, args, kwargs)
        self.pending_futures[task_id] = ft
        self.task_queue.append(task_tuple)
        return ft
    
    def submit(self, fn, *args, **kwargs):
        ft = self._make_append_task(fn, args, kwargs)
        # wake up the dispatch loop
        self._signal_thread(self._dispatch_thread_ctl_endpoint, socket_type=zmq.PUSH)
        return ft
    
    def submit_many(self, tasks):
        futures = [self._make_append_task(fn, args, kwargs) for (fn,args,kwargs) in tasks]
        self._signal_thread(self._dispatch_thread_ctl_endpoint, socket_type=zmq.PUSH)
        return futures
        
    def shutdown(self):
        if not self._shutdown_signaled:
            #self.master_announce_socket.send('shutdown')
            self._shutdown_signaled = True
            
            for endpoint in (self._dispatch_thread_ctl_endpoint,self._receive_thread_ctl_endpoint,self._announce_endpoint):
                self._signal_thread(endpoint, 'shutdown', socket_type=zmq.PUSH)            

class ZMQWMClient(ZMQBase):
    '''A multiplexing client, responsible for spawning worker processes,
    monitoring them, and forwarding traffic between workers and the server.
    Also will shut down if the server disappears.'''
    
    def __init__(self, n_workers,
                 upstream_task_endpoint, upstream_result_endpoint, upstream_announce_endpoint,
                 downstream_mode = 'ipc'):
        super(ZMQWMClient, self).__init__()
        
        self.n_workers = n_workers or multiprocessing.cpu_count()
        
        # cannot initialize context until after forks
        #self.context = zmq.Context.instance()
        self.context = None
        
        # Where we receive work, dispatch results, and receive announcements
        self.upstream_announce_endpoint = upstream_announce_endpoint
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_result_endpoint = upstream_result_endpoint

        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._announce_ctl_endpoint = 'inproc://_announce_{:x}'.format(id(self))
                
        # Queues of tasks and results
        #self.task_queue = deque()
        #self.pending_results = dict()
        
        # Information about worker processes
        self.workers = []
        
    def task_streamer_loop(self):
        '''Receive tasks from the server on behalf of clients.'''
        
    def result_streamer_loop(self):
        '''Dispatch results to the server on behalf of clients.'''
        
    def monitor_announcements_loop(self):
        '''Listen for announcements from the server (and internally about
        worker processes).'''
        
        ann_socket = self.context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')
        ann_socket.connect(self.upstream_announce_endpoint)
        
        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(self._announce_ctl_endpoint)
        
        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(ann_socket, zmq.POLLIN)
        
        last_announce = None
        remaining_interval = None
        try:
            while True:
                poll_results = dict(poller.poll(remaining_interval))
                if poll_results.get(ann_socket) == zmq.POLLIN:
                    message = ann_socket.recv()
                    if message == 'shutdown':
                        self.shutdown()
                        return
                    elif message == 'ping':
                        last_announce = time.time()
                        remaining_interval = self.server_heartbeat_interval * 1000
                                           
                now = time.time() 
                if now - last_announce < self.server_heartbeat_interval:
                    remaining_interval = (now - last_announce) * 1000
                else:
                    remaining_interval = self.server_heartbeat_interval * 1000
                
        finally:
            poller.unregister(ctlsocket)
            poller.unregister(ann_socket)
            ann_socket.close(linger=0)
            ctlsocket.close()
            
        
    def monitor_worker_loop(self):
        '''Spawn a worker process and monitor it.'''
        
    def startup(self):
        '''Fork worker processes and begin monitoring the server.'''
    
    def shutdown(self):
        '''Terminate worker processes and shut down'''
                

class ZMQWMClientWorker(ZMQBase):
    '''A worker process'''
    
    def __init__(self, upstream_task_endpoint, upstream_result_endpoint):
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_result_endpoint = upstream_result_endpoint
        
    def startup(self):
        # Become a session leader so that we can kill child processes with this
        # process by killing the process group
        os.setsid()
        os.setpgrp()
        self.context = zmq.Context()
        rc = 4
        try:
            rc = self._rpd_loop()
        finally:
            self.context.destroy(linger=0)
            del self.context
            sys.exit(rc)
            
    def _rpd_loop(self):
        '''Run the receive-perform-dispatch loop'''
        # reinitialize random number generator in this process
        random.seed()
        
        task_socket = self.context.socket(zmq.PULL)
        result_socket = self.context.socket(zmq.PUSH)
        task_socket.connect(self.upstream_task_endpoint)
        result_socket.connect(self.upstream_result_endpoint)
        associated_server = None
        try:
            while True:
                (server_id, task_id, fn, args, kwargs) = task_socket.recv_pyobj()
                if associated_server is None:
                    associated_server = server_id
                elif server_id != associated_server:
                    log.error('received task from {} when expecting task from {} (zombie process?); exiting'
                              .format(server_id, associated_server))
                    return 2
    
                try:
                    result = fn(*args, **kwargs)
                except Exception as e:
                    result_tuple = (server_id, task_id, 'exception', (e, traceback.format_exc()))
                else:
                    result_tuple = (server_id, task_id, 'result', result)
                result_socket.send_pyobj(result_tuple)
        except KeyboardInterrupt:
            log.info('received SIGINT; shutting down cleanly')
            return 0
        finally:
            task_socket.close(linger=0)
            result_socket.close(linger=0)


class ZMQWorkManager(ZMQWMServer,WorkManager):
    pass

        
atexit.register(ZMQBase.remove_ipc_endpoints)
      
