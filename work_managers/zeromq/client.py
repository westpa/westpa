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
Client for ZeroMQ work manager.

This coordinates the management of a number of processes, which use
ZeroMQ to communicate with the server.
'''

from __future__ import division, print_function; __metaclass__ = type

import os, sys, time, threading, multiprocessing, signal, traceback, logging, warnings
import zmq

import work_managers

# Import splat, so that we don't have to enumerate the constants we care
# about. Put the log instantiation afterwards, or else we will import log
# from core.
from core import *

from . import recvall, unpickle_frame, WorkerTerminated
from . import (DEFAULT_HANGCHECK_INTERVAL, DEFAULT_MAX_TASKQUEUE_SIZE, DEFAULT_SERVER_HEARTBEAT_INTERVAL, 
               DEFAULT_SHUTDOWN_TIMEOUT, DEFAULT_TASK_TIMEOUT, DEFAULT_TASKQUEUE_WAIT)

log = logging.getLogger(__name__)

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

    #NOTE: UNUSED (instead, see work_manager.py's add_args class method)
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
                                      assume the server has crashed and shut down. (Default: {} seconds.)'''
                                      .format(DEFAULT_SERVER_HEARTBEAT_INTERVAL))
        wm_group.add_argument(wmenv.arg_flag('zmq_task_timeout'), metavar='TIMEOUT', type=int,
                              help='''Kill worker processes that take longer than TIMEOUT seconds
                                      (default: no limit).''')
        wm_group.add_argument(wmenv.arg_flag('zmq_client_comm_mode'), metavar='MODE',
                              help='''Use the given MODE ('ipc' for Unix sockets and 'tcp' for
                              TCP/IP sockets) to communicate with worker processes. 'ipc' is more
                              efficient, but may not work properly if a node-local temporary filesystem
                              (e.g. /tmp) is not available. (Default: 'ipc')''')
                                             
    
    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 
        
        n_workers = wmenv.get_val('n_workers', multiprocessing.cpu_count(), int)
        task_timeout = wmenv.get_val('zmq_task_timeout', DEFAULT_TASK_TIMEOUT, int)
        heartbeat_interval = wmenv.get_val('zmq_heartbeat_interval', DEFAULT_SERVER_HEARTBEAT_INTERVAL, int)
        shutdown_timeout = wmenv.get_val('zmq_worker_shutdown_timeout', DEFAULT_SHUTDOWN_TIMEOUT, int)
        hangcheck_interval = wmenv.get_val('zmq_hangcheck_interval', DEFAULT_HANGCHECK_INTERVAL, int)
        client_comm_mode = wmenv.get_val('zmq_client_comm_mode')

        # if individual endpoints are named, we use these
        tests_old = [not bool(wmenv.get_val('zmq_task_endpoint')),
                 not bool(wmenv.get_val('zmq_result_endpoint')),
                 not bool(wmenv.get_val('zmq_announce_endpoint')),
                 not bool(wmenv.get_val('zmq_listen_endpoint'))]
        tests_new = [not bool(wmenv.get_val('zmq_upstream_task_endpoint')),
                 not bool(wmenv.get_val('zmq_upstream_result_endpoint')),
                 not bool(wmenv.get_val('zmq_upstream_announce_endpoint')),
                 not bool(wmenv.get_val('zmq_upstream_listen_endpoint'))]

        if all(tests_old) and all(tests_new):
            # No endpoints specified; use server/router info file
            upstream_info_filename = wmenv.get_val('zmq_read_info') or wmenv.get_val('zmq_info')
            if upstream_info_filename is None:
                raise EnvironmentError('neither endpoints nor upstream (server or router) info file specified')
            else:
                import json
                try:
                    upstream_info = json.load(open(upstream_info_filename,'rt'))
                    task_endpoint = upstream_info['task_endpoint']
                    result_endpoint = upstream_info['result_endpoint']
                    announce_endpoint = upstream_info['announce_endpoint'] 
                    listen_endpoint = upstream_info['listen_endpoint']   
                    
                except Exception as e:
                    raise EnvironmentError('cannot load upstream info file {!r}: {}'.format(upstream_info_filename,e))

        elif (not all(tests_old) and any(tests_old)) or (not all(tests_new) and any(tests_new)):
            raise ValueError('either none or all three server endpoints must be specified')
        #Use new-style, unambiguous endpoint args
        elif not all(tests_new):
            task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_upstream_task_endpoint'),allow_wildcard_host=False)
            result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_upstream_result_endpoint'),allow_wildcard_host=False)
            announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_upstream_announce_endpoint'),allow_wildcard_host=False)
            listen_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_upstream_listen_endpoint'),allow_wildcard_host=False)
        else:
            log.debug('using old style endpoint for client')
            task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_task_endpoint'),allow_wildcard_host=False)
            result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_result_endpoint'),allow_wildcard_host=False)
            announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_announce_endpoint'),allow_wildcard_host=False)
            listen_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_listen_endpoint'),allow_wildcard_host=False)
        
        return cls(task_endpoint, result_endpoint, announce_endpoint, listen_endpoint, server_heartbeat_interval = heartbeat_interval,
                   comm_mode = client_comm_mode, n_workers = n_workers, task_timeout = task_timeout,
                   shutdown_timeout = shutdown_timeout, hangcheck_interval = hangcheck_interval)
                     
    @classmethod
    def make_tcp_endpoint(cls):
        return 'tcp://127.0.0.1:{}'.format(randport())
    
    def __init__(self, upstream_task_endpoint, upstream_result_endpoint, upstream_announce_endpoint, upstream_update_endpoint,
                 server_heartbeat_interval = None, comm_mode = None, n_workers = None, task_timeout = None,
                 shutdown_timeout = None, hangcheck_interval = None):
        
        # Argument normalization
        server_heartbeat_interval = server_heartbeat_interval or DEFAULT_SERVER_HEARTBEAT_INTERVAL
        task_timeout = task_timeout or DEFAULT_TASK_TIMEOUT
        shutdown_timeout = shutdown_timeout or DEFAULT_SHUTDOWN_TIMEOUT
        if n_workers is None:
            n_workers = n_workers or multiprocessing.cpu_count()
        comm_mode = comm_mode or 'ipc'
        if not hangcheck_interval:
            hangcheck_interval = min(server_heartbeat_interval, task_timeout or sys.maxint)
        
        super(ZMQClient,self).__init__(server_heartbeat_interval)
        
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_result_endpoint = upstream_result_endpoint
        self.upstream_announce_endpoint = upstream_announce_endpoint
        self.upstream_listen_endpoint = upstream_update_endpoint
        
        self.context = None # this really shouldn't be instantiated until after forks, just to be safe
        
        self.n_workers = n_workers

        # How long we wait for worker processes to shutdown on SIGTERM
        # before moving on with SIGKILL
        self.shutdown_timeout = shutdown_timeout
        
        # How long we wait for a response from a worker process before
        # declaring it hung and terminating it -- None means do not check.
        self.worker_task_timeout = task_timeout
        
        # How often (in s) we check for a hung worker
        self.hangcheck_interval = hangcheck_interval
        
        if comm_mode == 'ipc':
            self.worker_task_endpoint = self.make_ipc_endpoint()
            self.worker_result_endpoint = self.make_ipc_endpoint()
        elif comm_mode == 'tcp':
            self.worker_task_endpoint = self.make_tcp_endpoint()
            self.worker_result_endpoint = self.make_tcp_endpoint()
                
        self.worker_lock = threading.RLock()
        self.workers = {} # mapping of PID to Process objects
        self.worker_active_tasks = {} # mapping of PID to Task objects (without payloads)
        self.worker_last_dispatch = {} # mapping of PID to time when last dispatch occurred
        
        self.associated_server_id = None
        
        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._taskfwd_ctl_endpoint = 'inproc://_taskfwd_ctl_{:x}'.format(id(self))
        self._rslfwd_ctl_endpoint  = 'inproc://_rslfwd_ctl_{:x}'.format(id(self))
        self._monitor_ctl_endpoint = 'inproc://_monitor_ctl_{:x}'.format(id(self)) #socket that gets server announcements
        self._update_ctl_endpoint = 'inproc://_update_ctl_{:x}'.format(id(self)) #socket that sends updates to server listen socket
        
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

                self._update_thread = threading.Thread(target=self._update_server_loop)
                self._update_thread.start()
                
                # Wait on all four threads starting up before continuing
                ctlsocket.recv()
                ctlsocket.recv()
                ctlsocket.recv()
                ctlsocket.recv()
    
            finally:
                ctlsocket.close()
            
    def _shutdown(self):
        if not self._shutdown_signaled:
            self._shutdown_signaled = True
            for endpoint in (self._monitor_ctl_endpoint, self._rslfwd_ctl_endpoint,
                             self._taskfwd_ctl_endpoint, self._update_ctl_endpoint):
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
        self._update_thread.join()
        
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
                waittime = min(self.hangcheck_interval,self.server_heartbeat_interval)*1000
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
            ctlsocket.close()
            log.debug('client: exiting monitor loop')

    def _update_server_loop(self):
        '''(complement to server's _listen_loop)
        Contains socket that updates the server on this client's status, including its number of active workers'''

        ctlsocket = self._make_signal_socket(self._update_ctl_endpoint) #Signals this thread to update server or shutdown
        upstream_update_socket = self.context.socket(zmq.PUSH)
        upstream_update_socket.connect(self.upstream_listen_endpoint)
        
        self._signal_thread(self._startup_ctl_endpoint) #signal successful startup

        #Not necessary; first time check in loop will fail (because last_update initially set to 0)
        self._signal_thread(self._update_ctl_endpoint, 'update') #signal self to send an intial update to server

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)

        last_update = 0
        remaining_interval = self.server_heartbeat_interval #time interval for sending server updates

        this_node_id = self.instance_id

        try:
            while True:
                poll_results = dict(poller.poll(remaining_interval*1000))

                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)

                    if MSG_SHUTDOWN in messages:
                        #alert server of shutdown and exit loop.
                        #note that this will send shutdown update to server, even if shutdown is *because*
                        #server shut down - an unneeded messge, but shouldn't matter.
                        upstream_update_socket.send_pyobj((MSG_SHUTDOWN, None, this_node_id, None)) 
                        return

                    elif 'update' in messages: #Update server

                        inner_pid_message = [] #worker pid's
                        inner_id_message = [] #worker id's

                        for pid, worker in self.workers.items():
                            inner_pid_message.append(pid)
                            inner_id_message.append(worker.instance_id)

                        inner_pid_message = tuple(inner_pid_message)
                        inner_id_message = tuple(inner_id_message)

                        assert len(inner_pid_message) == len(inner_id_message)

                        current_n_workers = len(inner_pid_message)
                        outer_message = (MSG_PING, None, this_node_id, current_n_workers)

                        upstream_update_socket.send_pyobj(outer_message, flags=zmq.SNDMORE)
                        upstream_update_socket.send_pyobj(inner_pid_message, flags=zmq.SNDMORE)
                        upstream_update_socket.send_pyobj(inner_id_message)

                    else:
                        for message in messages:
                            upstream_update_socket.send(message)

                now = time.time()
                if now - last_update >= self.server_heartbeat_interval:

                    last_update = now

                    self._signal_thread(self._update_ctl_endpoint, 'update') #signal it's time to update server



        finally:
            poller.unregister(ctlsocket)
            upstream_update_socket.close(linger=0)
            ctlsocket.close()
            log.debug('client: exiting update loop')
            
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
            ctlsocket.close()
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
            ctlsocket.close()
            log.debug('exiting result forwarding loop')

    @property
    def is_master(self):
        '''True if this is the master process for task distribution. This is necessary, e.g., for
        MPI, where all processes start identically and then must branch depending on rank.'''
        return False
