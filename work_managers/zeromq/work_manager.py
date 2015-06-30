'''
Created on Jun 10, 2015

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

from core import ZMQCore, Message, Task, Result, ZMQWorkerMissing, ZMQWMEnvironmentError, IsNode
from core import randport
from worker import ZMQWorker
from node import ZMQNode
import work_managers
from work_managers import WorkManager, WMFuture
import multiprocessing

from core import PassiveMultiTimer

import zmq

from collections import deque

import socket, re, json

class ZMQWorkManager(ZMQCore,WorkManager,IsNode):
    
    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env
            
        wm_group = parser.add_argument_group('options for ZeroMQ ("zmq") work manager (master or node)')
        
        wm_group.add_argument(wmenv.arg_flag('zmq_mode'), metavar='MODE', choices=('master','node','server','client'),
                              help='Operate as a master (server) or a node (workers/client). '
                                  +'"server" is a deprecated synonym for "master" and "client" is a '
                                  +'deprecated synonym for "node".')
        wm_group.add_argument(wmenv.arg_flag('zmq_comm_mode'), metavar='COMM_MODE', choices=('ipc', 'tcp'),
                              help='Use the given communication mode -- TCP or IPC (Unix-domain) -- sockets '
                                  +'for communication within a node. IPC (the default) may be more '
                                  +'efficient but is not available on (exceptionally rare) systems '
                                  +'without node-local storage (e.g. /tmp); on such systems, TCP may be used instead.')
        wm_group.add_argument(wmenv.arg_flag('zmq_write_host_info'), metavar='INFO_FILE',
                              help='Store hostname and port information needed to connect to this instance '
                                  +'in INFO_FILE. This allows the master and nodes assisting in '
                                  +'coordinating the communication of other nodes to choose ports '
                                  +'randomly. Downstream nodes read this file with '
                                  +wmenv.arg_flag('zmq_read_host_info') + ' and know where how to connect.')
        wm_group.add_argument(wmenv.arg_flag('zmq_read_host_info'), metavar='INFO_FILE',
                              help='Read hostname and port information needed to connect to the master '
                                  +'(or other coordinating node) from INFO_FILE. '
                                  +'This allows the master and nodes assisting in '
                                  +'coordinating the communication of other nodes to choose ports '
                                  +'randomly, writing that information with '
                                  +wmenv.arg_flag('zmq_write_host_info') + ' for this instance to read.')
        wm_group.add_argument(wmenv.arg_flag('zmq_upstream_rr_endpoint'), metavar='ENDPOINT',
                              help='ZeroMQ endpoint to which to send request/response (task and result) '
                                  +'traffic toward the master.')
        wm_group.add_argument(wmenv.arg_flag('zmq_upstream_ann_endpoint'), metavar='ENDPOINT',
                              help='ZeroMQ endpoint on which to receive announcement '
                                  +'(heartbeat and shutdown notification) traffic from the master.')
        wm_group.add_argument(wmenv.arg_flag('zmq_downstream_rr_endpoint'), metavar='ENDPOINT',
                              help='ZeroMQ endpoint on which to listen for request/response '
                                  +'(task and result) traffic from subsidiary workers.')
        wm_group.add_argument(wmenv.arg_flag('zmq_downstream_ann_endpoint'), metavar='ENDPOINT',
                              help='ZeroMQ endpoint on which to send announcement '
                                  +'(heartbeat and shutdown notification) traffic toward workers.')
        wm_group.add_argument(wmenv.arg_flag('zmq_master_heartbeat'), metavar='MASTER_HEARTBEAT',
                              type=float,
                              help='Every MASTER_HEARTBEAT seconds, the master announces its presence '
                                  +'to workers.')
        wm_group.add_argument(wmenv.arg_flag('zmq_worker_heartbeat'), metavar='WORKER_HEARTBEAT',
                              type=float,
                              help='Every WORKER_HEARTBEAT seconds, workers announce their presence '
                                  +'to the master.')
        wm_group.add_argument(wmenv.arg_flag('zmq_timeout_factor'), metavar='FACTOR',
                              type=float,
                              help='Scaling factor for heartbeat timeouts. '
                                  +"If the master doesn't hear from a worker in WORKER_HEARTBEAT*FACTOR, "
                                  +"the worker is assumed to have crashed. If a worker doesn't hear from "
                                  +"the master in MASTER_HEARTBEAT*FACTOR seconds, the master is assumed "
                                  +"to have crashed. Both cases result in shutdown. "
                                  )
        wm_group.add_argument(wmenv.arg_flag('zmq_startup_timeout'), metavar='STARTUP_TIMEOUT', 
                              type=float,
                              help='Amount of time (in seconds) to wait for communication between '
                                  +'the master and at least one worker. This may need to be changed '
                                  +'on very large, heavily-loaded computer systems that start all processes '
                                  +'simultaneously. '
                                  )
        wm_group.add_argument(wmenv.arg_flag('zmq_shutdown_timeout'), metavar='SHUTDOWN_TIMEOUT', 
                              type=float,
                              help='Amount of time (in seconds) to wait for workers to shut down.')
    
    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env
            
        # determine mode
        mode = wmenv.get_val('zmq_mode', 'master').lower()
        if mode in {'master', 'server'}:
            mode = 'master'
        elif mode in {'node', 'client'}:
            mode = 'node'
        else:
            raise ValueError('invalid ZMQ work manager mode {!r}'.format(mode))
        
        # determine number of workers
        # 0 with mode=='master' is a dedicated master
        # 0 with mode=='node' is a dedicated communications process (former ZMQRouter)
        n_workers = wmenv.get_val('n_workers', multiprocessing.cpu_count(), int)
        
        # We set this at the class level, because outside of testing, a node either
        # can support IPC or it can't, and there is no obvious need (currently)
        # to support both modes on an instance-by-instance basis
        comm_mode = wmenv.get_val('zmq_comm_mode', cls.default_comm_mode)
        ZMQWorkManager.internal_transport = comm_mode
        ZMQWorker.internal_transport = comm_mode
        ZMQNode.internal_transport = comm_mode
        
        write_host_info = wmenv.get_val('zmq_write_host_info')
        read_host_info  = wmenv.get_val('zmq_read_host_info')
        master_heartbeat = wmenv.get_val('zmq_master_heartbeat', cls.default_master_heartbeat, float)
        worker_heartbeat = wmenv.get_val('zmq_worker_heartbeat', cls.default_worker_heartbeat, float)
        timeout_factor = wmenv.get_val('zmq_timeout_factor', cls.default_timeout_factor, float)
        startup_timeout = wmenv.get_val('zmq_startup_timeout', cls.default_startup_timeout, float)
        
        
        if mode == 'master':
            instance = ZMQWorkManager(n_workers)
        else: # mode =='node'
            
            upstream_info = {}
            if read_host_info:
                upstream_info.update(cls.read_host_info(read_host_info))
            log.debug('upstream_info: {!r}'.format(upstream_info))
            
            upstream_rr_endpoint = wmenv.get_val('zmq_upstream_rr_endpoint', upstream_info.get('rr_endpoint'))
            upstream_ann_endpoint = wmenv.get_val('zmq_upstream_ann_endpoint', upstream_info.get('ann_endpoint'))
            
            if not (upstream_rr_endpoint and upstream_ann_endpoint):
                raise ZMQWMEnvironmentError('at least one upstream endpoint unspecified')
            
            # expand hostnames, if present, to IP addresses
            # reject wildcard hostnames, which is a logic error (can't connect to a host
            # without specifying an address)
            upstream_rr_endpoint = cls.canonicalize_endpoint(upstream_rr_endpoint, allow_wildcard_host=False)
            upstream_ann_endpoint = cls.canonicalize_endpoint(upstream_ann_endpoint, allow_wildcard_host=False)
            
            log.debug('upstream_rr_endpoint = {}'.format(upstream_rr_endpoint))
            log.debug('upstream_ann_endpoint = {}'.format(upstream_ann_endpoint))
            
            instance = ZMQNode(upstream_ann_endpoint=upstream_ann_endpoint, 
                               upstream_rr_endpoint=upstream_rr_endpoint, 
                               n_local_workers=n_workers)
        
        # Both server and node bind downstream endpoints, so that users get fan-out communications
        # "for free" when starting up a computational node    
        downstream_rr_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_rr_endpoint', 
                                                                         'tcp://*:{}'.format(randport())))
        downstream_ann_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_ann_endpoint', 
                                                                          'tcp://*:{}'.format(randport())))
        instance.downstream_rr_endpoint = downstream_rr_endpoint
        instance.downstream_ann_endpoint = downstream_ann_endpoint

        
        instance.master_beacon_period = master_heartbeat
        instance.worker_beacon_period = worker_heartbeat
        instance.timeout_factor = timeout_factor
        instance.startup_timeout = startup_timeout
        
        assert isinstance(instance, IsNode)
        for worker in instance.local_workers:
            worker.master_beacon_period = master_heartbeat
            worker.worker_beacon_period = worker_heartbeat
            worker.timeout_factor = timeout_factor
            worker.startup_timeout = startup_timeout
        
        # We always write host info (since we are always either master or node)
        # we choose not to in the special case that read_host_info is '' but not None
        # (None implies nothing found on command line or in environment variables, but ''
        # implies that it was found somewhere but it is empty)
        if write_host_info is not None and write_host_info != '':
            instance.write_host_info(write_host_info)

        log.debug('prepared {!r} with:'.format(instance))
        log.debug('n_workers = {}'.format(n_workers))
        for attr in ('master_beacon_period', 'worker_beacon_period', 'startup_timeout', 'timeout_factor',
                     'downstream_rr_endpoint', 'downstream_ann_endpoint'):
            log.debug('{} = {!r}'.format(attr, getattr(instance, attr)))
                
        return instance
    
    @classmethod
    def read_host_info(cls, filename):
        return json.load(open(filename,'rt'))
        
    

    @classmethod
    def canonicalize_endpoint(cls, endpoint, allow_wildcard_host = True):
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
    
    def __init__(self, n_local_workers=1):
        ZMQCore.__init__(self)
        WorkManager.__init__(self)
        IsNode.__init__(self, n_local_workers)
                                
        # Futures indexed by task ID
        self.futures = dict()
        
        # Tasks pending distribution
        self.outgoing_tasks = deque()
        
        # Tasks being processed by workers (indexed by worker_id)
        self.assigned_tasks = dict()
        
        # Identity information and last contact from workers
        self.worker_information = dict() # indexed by worker_id
        self.worker_timeouts = PassiveMultiTimer() # indexed by worker_id
                
        # Number of seconds between checks to see which workers have timed out
        self.worker_timeout_check = 5.0
        
        # Amount of time to wait for stray requests to arrive so that workers shut down properly
        self.shutdown_timeout = 0.5
        
        self.master_id = self.node_id
            
    @property
    def n_workers(self):
        return len(self.worker_information)
                
    def submit(self, fn, args=None, kwargs=None):
        if self.futures is None:
            # We are shutting down
            raise ZMQWMEnvironmentError('work manager is shutting down')
        future = WMFuture()
        task = Task(fn, args or (), kwargs or {}, task_id = future.task_id)
        self.futures[task.task_id] = future
        self.outgoing_tasks.append(task)
        # Wake up the communications loop (if necessary) to announce new tasks
        self.send_inproc_message(Message.TASKS_AVAILABLE)
        return future

    def submit_many(self, tasks):
        if self.futures is None:
            # We are shutting down
            raise ZMQWMEnvironmentError('work manager is shutting down')
        futures = []        
        for (fn,args,kwargs) in tasks:
            future = WMFuture()
            task = Task(fn, args, kwargs, task_id = future.task_id)
            self.futures[task.task_id] = future
            self.outgoing_tasks.append(task)
            futures.append(future)
        # Wake up the communications loop (if necessary) to announce new tasks            
        self.send_inproc_message(Message.TASKS_AVAILABLE)
        return futures

    def send_message(self, socket, message, payload=None, flags=0):
        message = Message(message, payload)
        message.master_id = self.node_id
        super(ZMQWorkManager,self).send_message(socket, message, payload, flags)
        
    def handle_result(self, socket, msg):
        self.send_ack(socket,msg)
        with self.message_validation(msg):
            assert msg.message == Message.RESULT
            assert isinstance(msg.payload, Result)
            assert msg.payload.task_id in self.futures
            assert self.assigned_tasks[msg.src_id].task_id == msg.payload.task_id
                        
        result = msg.payload
        
        future = self.futures.pop(result.task_id)
        del self.assigned_tasks[msg.src_id]
        if result.exception is not None:
            future._set_exception(result.exception, result.traceback)
        else:
            future._set_result(result.result)
            
    def handle_task_request(self, socket, msg):
        if not self.outgoing_tasks:
            # No tasks available
            self.send_nak(socket,msg)
        else:
            task = self.outgoing_tasks.popleft()
            
            worker_id = msg.src_id
            self.assigned_tasks[worker_id] = task
            
            self.send_message(socket, Message.TASK, task)
            
    def update_worker_information(self, msg):
        if msg.message == Message.IDENTIFY:
            with self.message_validation(msg):
                assert isinstance(msg.payload, dict)
            self.worker_information[msg.src_id] = msg.payload
        else:
            self.worker_information[msg.src_id] = {}
        
        try:
            self.worker_timeouts.reset(msg.src_id)
        except KeyError:
            self.worker_timeouts.add_timer(msg.src_id,self.worker_beacon_period*self.timeout_factor)
        
    def check_workers(self):
        expired_worker_ids = self.worker_timeouts.which_expired()
        for expired_worker_id in expired_worker_ids:
            try:
                worker_description = '{!s} ({!s})'.format(expired_worker_id, 
                                                          self.worker_information[expired_worker_id]['description'])
            except KeyError:
                worker_description = str(expired_worker_id)
            
            self.log.error('no contact from worker {}'.format(expired_worker_id, worker_description))
               
            self.remove_worker(expired_worker_id)            
                                        
    def remove_worker(self, worker_id):
        try:
            expired_task = self.assigned_tasks.pop(worker_id)
        except KeyError:
            pass
        else:
            self.log.error('aborting task {!r} running on expired worker {!s}'
                           .format(expired_task, worker_id))
            future = self.futures.pop(expired_task.task_id)
            future._set_exception(ZMQWorkerMissing('worker running this task disappeared'))
        del self.worker_information[worker_id]
        
    def shutdown_clear_tasks(self):
        '''Abort pending tasks with error on shutdown.'''
        while self.futures:
            task_id, future = self.futures.popitem()
            future._set_exception(ZMQWMEnvironmentError('work manager shut down during task'))
        self.futures = None
            
    
    def comm_loop(self):
        self.context = zmq.Context()
        
        rr_socket = self.context.socket(zmq.REP)
        ann_socket = self.context.socket(zmq.PUB)
        
        for endpoint in (self.local_rr_endpoint, self.downstream_rr_endpoint):
            if endpoint: rr_socket.bind(endpoint)
            
        for endpoint in (self.local_ann_endpoint, self.downstream_ann_endpoint):
            if endpoint: ann_socket.bind(endpoint)

        inproc_socket = self.context.socket(zmq.SUB)
        inproc_socket.setsockopt(zmq.SUBSCRIBE,'')
        inproc_socket.bind(self.inproc_endpoint)
        
        poller = zmq.Poller()
        poller.register(inproc_socket, zmq.POLLIN)
        poller.register(rr_socket, zmq.POLLIN)
        
        timers = PassiveMultiTimer()
        timers.add_timer('tasks_avail', self.master_beacon_period)
        timers.add_timer('master_beacon', self.master_beacon_period)
        timers.add_timer('worker_timeout_check', self.worker_beacon_period*self.timeout_factor)
        timers.add_timer('startup_timeout', self.startup_timeout)
        timers.reset()
        
        self.log.debug('master beacon period: {!r}'.format(self.master_beacon_period))
        self.log.debug('startup timeout: {!r}'.format(self.startup_timeout))
                
        peer_found = False
        
        try:
            # Send a master alive message immediately; it will get discarded if necessary
            self.send_message(ann_socket, Message.MASTER_BEACON)
            
            while True:
                # If a timer is already expired, next_expiration_in() will return 0, which
                # zeromq interprets as infinite wait; so instead we select a 1 ms wait in this
                # case.
                timeout = (timers.next_expiration_in() or 0.001)*1000
                # Wake up every second to check for signals
                timeout = min(timeout, 1000)
                poll_results = dict(poller.poll(timeout))
                                
                if inproc_socket in poll_results:
                    msgs = self.recv_all(inproc_socket,validate=False)
                    # Check for shutdown; do nothing else if shutdown is signalled
                    if Message.SHUTDOWN in (msg.message for msg in msgs):
                        self.log.debug('shutdown received')
                        break
                    # Check for any other wake-up messages
                    for msg in msgs:
                        if msg.message == Message.TASKS_AVAILABLE:
                            self.send_message(ann_socket, Message.TASKS_AVAILABLE)
                
                if rr_socket in poll_results:
                    msg = self.recv_message(rr_socket)
                    self.update_worker_information(msg)
                    
                    if msg.message == Message.TASK_REQUEST:
                        self.handle_task_request(rr_socket, msg)
                    elif msg.message == Message.RESULT:
                        self.handle_result(rr_socket, msg)
                    else:
                        self.send_ack(rr_socket, msg)
                        
                    if self.worker_information:
                        peer_found = True
                
                if timers.expired('tasks_avail'):
                    if self.outgoing_tasks:
                        self.send_message(ann_socket, Message.TASKS_AVAILABLE)
                    timers.reset('tasks_avail')
                    
                if timers.expired('master_beacon'):
                    self.send_message(ann_socket, Message.MASTER_BEACON)
                    timers.reset('master_beacon')
                    
                if peer_found and timers.expired('worker_timeout_check'):
                    self.check_workers()
                    if not self.worker_information:
                        self.log.error('all workers disappeared; exiting')
                        break
                    timers.reset('worker_timeout_check')
                    
                if not peer_found and timers.expired('startup_timeout'):
                    self.log.error('startup phase elapsed with no contact from workers; shutting down')
                    while self.futures:
                        future = self.futures.popitem()[1]
                        future._set_exception(ZMQWorkerMissing('no workers available'))
                    break
                
            # Post a shutdown message
            self.log.debug('sending shutdown on ann_socket')
            self.send_message(ann_socket, Message.SHUTDOWN)            
            poller.unregister(inproc_socket)
            
            # Clear tasks
            self.shutdown_clear_tasks()
            
            # Clear incoming queue of requests, to let clients exit request/reply states gracefully
            # (clients will still timeout in these states if necessary)
            timers.add_timer('shutdown', self.shutdown_timeout)
            while not timers.expired('shutdown'):
                poll_results = dict(poller.poll(self.shutdown_timeout / 10 * 1000))
                if rr_socket in poll_results:
                    msg = self.recv_message(rr_socket)
                    self.send_nak(rr_socket, msg)

        finally:
            self.context.destroy(linger=1)
            self.context = None
            self.remove_ipc_endpoints()
    
    def startup(self):
        IsNode.startup(self)
        super(ZMQWorkManager,self).startup()
        
    def shutdown(self):
        self.signal_shutdown()
        IsNode.shutdown(self)
        self.join()
        super(ZMQWorkManager,self).shutdown()
        