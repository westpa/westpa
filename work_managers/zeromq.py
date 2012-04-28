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
    
class ZMQWorkManager(WorkManager):
    _ipc_endpoints = []
    
    def __init__(self, n_workers = None):
        super(ZMQWorkManager,self).__init__()
        
        # Number of worker processes
        self.n_workers = n_workers or multiprocessing.cpu_count()
                
        # ZeroMQ context
        self.context = None
        
        # number of seconds between announcements of where to connect to the master
        self.master_heartbeat_interval = 10
        
        # This hostname
        self.hostname = socket.gethostname()
        self.host_id = '{:s}-{:d}'.format(self.hostname, os.getpid())
        self.instance_id = uuid.uuid4()
                
        self._ipc_endpoints = []

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
        raise NotImplementedError('this work manager must be made a master or worker before startup')
    
    def shutdown(self):
        raise NotImplementedError
    
    def __enter__(self):
        self.startup()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_traceback):
        self.shutdown()
        return False
    
SHUTDOWN_SENTINEL = ('shutdown',)

class ZMQWMMaster(ZMQWorkManager):
    
    # tasks are tuples of (task_id, function, args, keyword args)
    # results are tuples of (task_id, {'result', 'exception'}, value)
    
    def __init__(self, n_workers, 
                 master_announce_endpoint, master_task_endpoint, master_result_endpoint):
        super(ZMQWMMaster,self).__init__(n_workers)
        
        self.context = zmq.Context.instance()
        
        # Where we send out announcements
        self.master_announce_endpoint = master_announce_endpoint
        
        # where we send out work
        self.master_task_endpoint = master_task_endpoint
        
        # Where we receive results
        self.master_result_endpoint = master_result_endpoint
 
                        
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
        
    def _signal_thread(self, endpoint, message='', socket_type=zmq.PUSH):
        socket = self.context.socket(socket_type)
        socket.connect(endpoint)
        socket.send(message)
        socket.close()
        del socket

                        
    def _dispatch_loop(self):
        # a 1-1 mapping between items added to the task queue and ZMQ sends
        
        # Bind the task distributor socket
        master_task_socket = self.context.socket(zmq.PUSH)
        master_task_socket.setsockopt(zmq.HWM,1)
        master_task_socket.bind(self.master_task_endpoint)

        # Create a control socket to wake up the loop        
        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(self._dispatch_thread_ctl_endpoint)
        
        #self._signal_startup_ctl()
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
        remaining_interval = self.master_heartbeat_interval
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
                if now - last_announce < self.master_heartbeat_interval:
                    remaining_interval = now - last_announce
                else:
                    remaining_interval = self.master_heartbeat_interval
                
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
        
atexit.register(ZMQWorkManager.remove_ipc_endpoints)
      
