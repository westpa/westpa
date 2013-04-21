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

from __future__ import division, print_function; __metaclass__ = type

import time, threading, logging
import Queue as queue

import zmq

from work_managers import WMFuture
from . import recvall, pickle_to_frame, unpickle_frame
from . import (DEFAULT_HANGCHECK_INTERVAL, DEFAULT_MAX_TASKQUEUE_SIZE, DEFAULT_SERVER_HEARTBEAT_INTERVAL, 
               DEFAULT_SHUTDOWN_TIMEOUT, DEFAULT_TASK_TIMEOUT, DEFAULT_TASKQUEUE_WAIT)

from core import *


log = logging.getLogger(__name__)

# A singleton object to put on the server task queue when shutting down,
# to break any blocking waits on it.
_ShutdownSentinel = object()

# Server-side task information
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

def _recvall_multipart(socket):
    #similar to 'recvall' package method, only works for multipart messages (for the listen socket)

    frames = []

    while True:
        try:
            frames.append(socket.recv_multipart(copy=False, flags=zmq.NOBLOCK))
        except zmq.ZMQError as err:
            if err.errno == zmq.EAGAIN:
                return frames
            else:
                raise               

class ZMQServer(ZMQBase):
    
    def __init__(self, master_task_endpoint, master_result_endpoint, master_announce_endpoint, master_listen_endpoint,
                 server_heartbeat_interval=DEFAULT_SERVER_HEARTBEAT_INTERVAL,
                 max_taskqueue_size=DEFAULT_MAX_TASKQUEUE_SIZE,
                 taskqueue_wait=DEFAULT_TASKQUEUE_WAIT):
        
        super(ZMQServer, self).__init__(server_heartbeat_interval)
        
        self.context = zmq.Context()

        self.clients = {} # Dictionary of client node ids => integer number of workers
                          # Sum of all values is the total number of workers responding to this server
                        
        # where we send out work
        self.master_task_endpoint = master_task_endpoint
        
        # Where we receive results
        self.master_result_endpoint = master_result_endpoint
 
        # Where we send out announcements
        self.master_announce_endpoint = master_announce_endpoint

        # Listen for updates from clients
        self.master_listen_endpoint = master_listen_endpoint
        
        # tasks awaiting dispatch
        self.task_queue = queue.Queue(max_taskqueue_size or 0)
        
        # upon request from a client, how long we wait for tasks to become
        # available before returning a "task unavailable" message 
        self.task_queue_wait = taskqueue_wait
        
        # tasks awaiting return, keyed by task ID
        self.pending_tasks = dict()
        
        self._shutdown_signaled = False

        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._dispatch_thread_ctl_endpoint = 'inproc://_dispatch_thread_ctl_{:x}'.format(id(self))        
        self._receive_thread_ctl_endpoint = 'inproc://_receive_thread_ctl_{:x}'.format(id(self))
        self._announce_endpoint = 'inproc://_announce_{:x}'.format(id(self))
        self._listen_endpoint = 'inproc://_listen_{:x}'.format(id(self))

    @property
    def n_workers(self):
        '''returns all the workers for all client processes'''
        return sum(self.clients.values())

    def startup(self):
        # start up server threads, blocking until their sockets are ready
        
        # create an inproc socket to sequence the startup of worker threads
        # each thread needs to write an empty message to this endpoint so
        # that startup() doesn't exit until all required sockets are open
        # and listening
        
        ctlsocket = self._make_signal_socket(self._startup_ctl_endpoint)
        
        #proper use here is to start a thread, then recv, and in the thread func
        #use _signal_startup_ctl() once all its sockets are bound
        self._dispatch_thread = threading.Thread(target=self._dispatch_loop)
        self._dispatch_thread.start()        
        
        self._receive_thread = threading.Thread(target=self._receive_loop)
        self._receive_thread.start()
        
        self._announce_thread = threading.Thread(target=self._announce_loop)
        self._announce_thread.start()
        
        self._listen_thread = threading.Thread(target=self._listen_loop)
        self._listen_thread.start()

        # These three recvs, together, block until all three threads
        # have started and the inproc communications endpoints are
        # bound.
        ctlsocket.recv() # dispatch
        ctlsocket.recv() # receive
        ctlsocket.recv() # announce
        ctlsocket.recv() # listen
        
        ctlsocket.close()
        
                        
    def _dispatch_loop(self):
        
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
                            task = self.task_queue.get(block=True, timeout=self.task_queue_wait)
                        except queue.Empty:
                            log.debug('server: no task available')
                            master_task_socket.send_pyobj((MSG_TASK_UNAVAILABLE, this_server_id, client_id, None))
                        else:
                            if task is _ShutdownSentinel:
                                return
                            
                            # Send task as a two-part message: a (small) metadata header
                            # and a (possibly large) payload
                            log.debug('server: dispatching task {!s} ({!r})'.format(task.task_id, task.fn))
                            master_task_socket.send_pyobj((MSG_TASK_AVAILABLE, this_server_id, client_id, task.task_id),
                                                          flags=zmq.SNDMORE)
                            master_task_socket.send(pickle_to_frame((task.fn, task.args, task.kwargs)), copy=False)
                            del task
                    else:
                        log.error('unknown/unsupported message received on task socket: {!r}'.format(tag))
                    del message, tag, server_id, client_id, _payload
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
                            task = self.pending_tasks.pop(task_id)
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
                                del retval
                            else:
                                log.error('unknown result type received for task {!s} ({!r})'.format(task_id, task.fn))
                            del result_type, result_payload, task
                    else:
                        log.error('unknown/unsupported message received on result socket: {!r}'.format(tag))
                    del frames, message, tag, server_id, client_id, payload
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

    def _listen_loop(self):
        #Thread for socket to 'listen' to downstream announcements from clients (for instance, how many workers each client has)
        #Kind of a reverse announce socket
        
        master_listen_socket = self.context.socket(zmq.PULL)
        master_listen_socket.bind(self.master_listen_endpoint)

        ctlsocket = self._make_signal_socket(self._listen_endpoint)

        self._signal_thread(self._startup_ctl_endpoint)

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(master_listen_socket, zmq.POLLIN)

        try:
            while True:
                poll_results = dict(poller.poll())

                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if MSG_SHUTDOWN in messages:
                        return

                if poll_results.get(master_listen_socket) == zmq.POLLIN:
                    messages = _recvall_multipart(master_listen_socket)

                    for frames in messages:
                        ##Unwrap the header of the message.  Note I'm not checking the rest (inner parts) of the message as of now
                        # (i.e. worker pid's and id's)

                        header = unpickle_frame(frames[0])

                        (tag, _, client_id, payload) = header

                        #Client node has shutdown
                        if tag == MSG_SHUTDOWN:
                            
                            try:
                                log.debug('sever: got shutdown message from client {!r}, removing client'.format(client_id))
                                self.clients.pop(client_id)
                            except KeyError:
                                log.error('server: received message from client that server doesn\'t know about - rogue client?')

                        #Clent node updating server
                        # may want to start a client-specific dict of last update times to decide to remove a client
                        # if it hasn't pinged server in specified time interval
                        elif tag == MSG_PING:

                            log.debug('server: got update from client {!r}'.format(client_id))
                            #new client started up
                            if client_id not in self.clients.keys():
                                log.debug('server: adding new client {!r} with {!r} workers'.format(client_id, payload))
                            else:
                                old_n_workers = self.clients[client_id]

                                log.debug('server: client {!r} update: used to have {!r} workers, now has {!r}'.format(client_id,
                                                                                                                       old_n_workers,
                                                                                                                       payload))
                            
                            self.clients[client_id] = payload

                        del header, frames

        finally:
            poller.unregister(ctlsocket)
            poller.unregister(master_listen_socket)
            master_listen_socket.close(linger=0)
            ctlsocket.close()
            log.debug('server: exiting listen loop')


            
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
            
            for endpoint in (self._dispatch_thread_ctl_endpoint,self._receive_thread_ctl_endpoint,
                             self._announce_endpoint, self._listen_endpoint):
                self._signal_thread(endpoint, MSG_SHUTDOWN)
                
            # Put a sentinel on the task queue to wake up waits on it
            self.task_queue.put(_ShutdownSentinel)
