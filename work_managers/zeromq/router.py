'''
A script containing a router to forward tasks/results between client(s) and server, and process announcements

Right now, implemented using three separate threads - not using devices for now
'''

import zmq, threading, time, logging, sys, socket, os, json, re, atexit
from . import pickle_to_frame, unpickle_frame, randport, recvall
from core import *
__metaclass__ = type

log = logging.getLogger(__name__)

from . import (DEFAULT_HANGCHECK_INTERVAL, DEFAULT_MAX_TASKQUEUE_SIZE, DEFAULT_SERVER_HEARTBEAT_INTERVAL, 
               DEFAULT_SHUTDOWN_TIMEOUT, DEFAULT_TASK_TIMEOUT, DEFAULT_TASKQUEUE_WAIT)

import work_managers

class ZMQDevice(ZMQBase, threading.Thread):
    '''
    homebrewed Device class (for forwarding result and task messages)

    Parameters
    ----------
    name - String name identifying this device
    upstream_endpoint - endpoint for upstream REQ socket
    downstream_endpoint - endpoint for downstream REP socket
    router_startup_ctl_endpoint - endpoint to signal router of succesful startup
    context - a zmq Context object
    server_heartbeat_interval - for initializing ZMQBase
    upstream_type - type of upstream socket (default is zmq.REQ)
    downstream_type - type of downstream socket (default is zmq.REP)
    '''

    def __init__(self, name, upstream_endpoint, downstream_endpoint, router_startup_ctl_endpoint, context = None,
                 upstream_type = zmq.REQ, downstream_type = zmq.REP):

        #set up arg
        context = context or zmq.Context() #Share context with router or no?

        ZMQBase.__init__(self)
        threading.Thread.__init__(self)

        self.setDaemon(1) # do I want this to be a daemon?

        self.name = name

        self.context = context 
        self.upstream_endpoint = upstream_endpoint
        self.downstream_endpoint = downstream_endpoint

        self._router_ctl_endpoint = router_startup_ctl_endpoint #endpoint to router object's ctlsocket

        self.upstream_type = upstream_type
        self.downstream_type = downstream_type

        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._shutdown_ctl_endpoint = 'inproc://_shutdown_ctl_{:x}'.format(id(self))

        self._shutdown_signalled = False
        self._running = False



    def __repr__(self):
        return "<Device {:d}: {!s}>".format(id(self), self.name)

    def __str__(self):
        return "Device '{!s}'".format(self.name)

    def startup(self):
        '''
        startup device thread
        '''
        
        ctlsocket = self._make_signal_socket(self._startup_ctl_endpoint) #startup ctlsocket

        self.start() # start up the device thread!

        ctlsocket.recv() #block until successful startup

        self._running = True
        self._shutdown_signalled = False
        log.debug('device: {} started up successful'.format(self.name))

        ctlsocket.close()

        


    def run(self): 
        #Contains main REQ/REP work
        #Forwards a message upstream, then returns the reply downstream
        #doesn't check anything (message headers, intended locations, etc.)

        ctlsocket = self._make_signal_socket(self._shutdown_ctl_endpoint) #A ctlsocket to receive a shutdown signal

        self._signal_thread(self._startup_ctl_endpoint) #signal successful startup

        #signal the router that startup is successful
        self._signal_thread(self._router_ctl_endpoint)

        upstream_socket = self.context.socket(self.upstream_type) #A (default REQ) socket to talk to upstream connection
        upstream_socket.connect(self.upstream_endpoint)

        downstream_socket = self.context.socket(self.downstream_type) #A (default REP) socket to talk to downstream connection
        downstream_socket.bind(self.downstream_endpoint)

        #Receive requests from downstream
        downstream_poller = zmq.Poller()
        downstream_poller.register(ctlsocket, zmq.POLLIN)
        downstream_poller.register(downstream_socket, zmq.POLLIN)

        #Receive replies from upstream
        upstream_poller = zmq.Poller()
        upstream_poller.register(ctlsocket, zmq.POLLIN)
        upstream_poller.register(upstream_socket, zmq.POLLIN)

        #Main work loop
        try:
            while True:

                #anything from downstream?
                downstream_poll_results = dict(downstream_poller.poll())

                if downstream_poll_results.get(ctlsocket) == zmq.POLLIN:
                    messages = recvall(ctlsocket)
                    if MSG_SHUTDOWN in messages:
                        return

                #Forward the message upstream
                if downstream_poll_results.get(downstream_socket) == zmq.POLLIN:
                    
                    request = downstream_socket.recv_multipart(copy=False) #Get downstream request
                    upstream_socket.send_multipart(request) #...and forward request upstream

                    #log.debug('device: sending request upstream')

                    upstream_poll_results = dict(upstream_poller.poll()) #Wait for upstream reply

                    if upstream_poll_results.get(ctlsocket) == zmq.POLLIN:
                        messages = recvall(ctlsocket)
                        if MSG_SHUTDOWN in messages:
                            return

                    #Forward reply downstream
                    if upstream_poll_results.get(upstream_socket) == zmq.POLLIN:

                        reply = upstream_socket.recv_multipart(copy=False) #Get reply from upstream
                        downstream_socket.send_multipart(reply)

                        #log.debug('device: sending reply downstream')

        #clean up
        finally:

            downstream_poller.unregister(ctlsocket)
            downstream_poller.unregister(downstream_socket)

            upstream_poller.unregister(ctlsocket)
            upstream_poller.unregister(upstream_socket)

            ctlsocket.close()
            upstream_socket.close()
            downstream_socket.close()

            log.debug('device: exiting run loop')

    def shutdown(self):
        '''shutdown this device''' 

        if not self._shutdown_signalled:

            self._shutdown_signalled = True
            self._running = False
            self._signal_thread(self._shutdown_ctl_endpoint, MSG_SHUTDOWN)

            log.debug("device: shutting down - bye")


class ZMQRouter(ZMQBase):


    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 

        wm_group = parser.add_argument_group('options for ZeroMQ ("zmq") router')
        wm_group.add_argument(wmenv.arg_flag('zmq_read_info'), metavar='UPSTREAM_INFO_FILE',
                              help='Read upstream endpoints from UPSTREAM_INFO_FILE for this router to connect to. '
                                   'This is helpful if running server and clients and routers on multiple '
                                   'machines which share a filesystem, as explicit hostnames/ports are not required')
        wm_group.add_argument(wmenv.arg_flag('zmq_write_info'), metavar='DOWNSTREAM_INFO_FILE',
                              help='Store this router\'s downstream endpoints in DOWNSTREAM_INFO_FILE '
                                   'for downstream clients to connect to. '
                                   'This is helpful if running server and clients and routers on multiple '
                                   'machines which share a filesystem, as explicit hostnames/ports are not required')
        wm_group.add_argument(wmenv.arg_flag('zmq_upstream_task_endpoint'), metavar='UPSTREAM_TASK_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint to connect upstream for task distribution''')
        wm_group.add_argument(wmenv.arg_flag('zmq_upstream_result_endpoint'), metavar='UPSTREAM_RESULT_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint to send results upstream''')
        wm_group.add_argument(wmenv.arg_flag('zmq_upstream_announce_endpoint'), metavar='UPSTREAM_ANNOUNCE_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint to connect upstream to receive announcements''')
        wm_group.add_argument(wmenv.arg_flag('zmq_downstream_task_endpoint'), metavar='DOWNSTREAM_TASK_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for downstream task distribution''')
        wm_group.add_argument(wmenv.arg_flag('zmq_downstream_result_endpoint'), metavar='DOWNSTREAM_RESULT_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint for downstream result collection''')
        wm_group.add_argument(wmenv.arg_flag('zmq_downstream_announce_endpoint'), metavar='DOWNSTREAM_ANNOUNCE_ENDPOINT',
                              help='''Use the given ZeroMQ endpoint to send announcements downstream''')
        wm_group.add_argument(wmenv.arg_flag('zmq_heartbeat_interval'), metavar='INTERVAL',
                              help='''If router has not
                                      heard from the server in approximately INTERVAL seconds, the router will
                                      assume the server has crashed and shut down itself and any downstream connections. 
                                      (Default: {} seconds.)'''.format(DEFAULT_SERVER_HEARTBEAT_INTERVAL))


    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 
        
        heartbeat_interval = wmenv.get_val('zmq_heartbeat_interval', DEFAULT_SERVER_HEARTBEAT_INTERVAL, int)

        ##SETUP UPSTREAM ENDPOINTS##
        #(UPSTREAM Endpoints for this router to connect to)
        tests = [not bool(wmenv.get_val('zmq_upstream_task_endpoint')),
                 not bool(wmenv.get_val('zmq_upstream_result_endpoint')),
                 not bool(wmenv.get_val('zmq_upstream_announce_endpoint'))]
        if all(tests):
            # No endpoints specified; use upstream/server info file
            upstream_info_filename = wmenv.get_val('zmq_read_info')
            if upstream_info_filename is None:
                raise EnvironmentError('neither upstream endpoints nor upstream server info file specified for router')
            else:
                
                try:
                    upstream_info = json.load(open(upstream_info_filename,'rt'))
                    upstream_task_endpoint = upstream_info['task_endpoint']
                    upstream_result_endpoint = upstream_info['result_endpoint']
                    upstream_announce_endpoint = upstream_info['announce_endpoint']    
                except Exception as e:
                    raise EnvironmentError('cannot load upstream info file {!r}: {}'.format(upstream_info_filename,e))                
        elif any(tests):
            raise ValueError('either none or all three upstream endpoints must be specified')
        else:
            upstream_task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_upstream_task_endpoint'),allow_wildcard_host=False)
            upstream_result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_upstream_result_endpoint'),allow_wildcard_host=False)
            upstream_announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_upstream_announce_endpoint'),allow_wildcard_host=False)

        ##SETUP DOWNSTREAM ENDPOINTS##
        #(Endpoints for DOWNSTREAM CLIENTS to connect to)
        router_info_filename = wmenv.get_val('zmq_write_info', 'zmq_router_info_{}.json'.format(uuid.uuid4().hex))
        
        # if individual endpoints are named, we use these
        tests = [not bool(wmenv.get_val('zmq_downstream_task_endpoint')),
                 not bool(wmenv.get_val('zmq_downstream_result_endpoint')),
                 not bool(wmenv.get_val('zmq_downstream_announce_endpoint'))]
        if all(tests):
            # Choose random ports
            downstream_task_endpoint = cls.canonicalize_endpoint('tcp://*')
            downstream_result_endpoint = cls.canonicalize_endpoint('tcp://*')
            downstream_announce_endpoint = cls.canonicalize_endpoint('tcp://*')
        elif any(tests):
            raise ValueError('either none or all three downstream endpoints must be specified')
        else:
            downstream_task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_task_endpoint'))
            downstream_result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_result_endpoint'))
            downstream_announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_announce_endpoint'))
        
        #Return ZMQRouter instance
        return cls(upstream_task_endpoint, upstream_result_endpoint, upstream_announce_endpoint, 
                   downstream_task_endpoint, downstream_result_endpoint, downstream_announce_endpoint,
                   server_heartbeat_interval=heartbeat_interval, router_info_filename=router_info_filename)

    def remove_router_info_file(self):
        if self.router_info_filename:
            filename = self.router_info_filename
            try:
                os.unlink(filename)
            except OSError as e:
                log.debug('could not remove router info file {!r}: {}'.format(filename, e))
            else:
                log.debug('removed router info file {!r}'.format(filename))

    def write_router_info(self, filename=None):
        '''Like server info files, this only writes the endpoints that this Router binds (i.e. endpoints to which downstream clients connect)'''

        filename = filename or self.router_info_filename
        hostname = socket.gethostname()
        with open(filename, 'wt') as infofile:
            log.debug('writing router info file {!r}'.format(filename))
            json.dump({'task_endpoint': re.sub(r'\*', hostname, self.downstream_task_endpoint),
                       'result_endpoint': re.sub(r'\*', hostname, self.downstream_result_endpoint),
                       'announce_endpoint': re.sub(r'\*', hostname, self.downstream_announce_endpoint)},
                      infofile)
        os.chmod(filename, 0600)


    def __init__(self, upstream_task_endpoint, upstream_result_endpoint, upstream_announce_endpoint, 
                 downstream_task_endpoint, downstream_result_endpoint, downstream_announce_endpoint,
                 server_heartbeat_interval=None, write_router_info=True, router_info_filename=None):

        

        #standardize heartbeat interval
        server_heartbeat_interval = server_heartbeat_interval or DEFAULT_SERVER_HEARTBEAT_INTERVAL


        super(ZMQRouter, self).__init__(server_heartbeat_interval)

        ##These should be initialized to the same values as they are for server and client

        self.context = zmq.Context()

        #Endpoints for server connections
        self.upstream_task_endpoint = upstream_task_endpoint 
        self.upstream_result_endpoint = upstream_result_endpoint 
        self.upstream_announce_endpoint = upstream_announce_endpoint 

        #Endpoints for client connections
        self.downstream_task_endpoint = downstream_task_endpoint
        self.downstream_result_endpoint = downstream_result_endpoint
        self.downstream_announce_endpoint = downstream_announce_endpoint

        ##Router info file##
        if write_router_info:
            self.router_info_filename = router_info_filename or 'zmq_write_info_{}.json'.format(uuid.uuid4().hex)
            self.write_router_info(self.router_info_filename)
            atexit.register(self.remove_router_info_file)
        else:
            self.router_info_filename = None

        self._shutdown_signalled = False
        self._running = False

        #ctlsocket endpoints for thread signalling
        self._startup_ctl_endpoint = 'inproc://_startup_ctl_{:x}'.format(id(self))
        self._announce_ctl_endpoint = 'inproc://_announce_ctl_{:x}'.format(id(self))

        #endpoint to receive device startup signals
        self._device_startup_ctl_endpoint = 'inproc://_device_startup_ctl_{:x}'.format(id(self))

        #Instance references to task and result devices
        self._task_device = None
        self._result_device = None

    def startup(self): 
        '''
        Starts up routing threads, blocks until all have started successfully
        '''

        ann_ctlsocket = self._make_signal_socket(self._startup_ctl_endpoint) #ctlsocket to receive successful announce thread startup signal
        device_ctlsocket = self._make_signal_socket(self._device_startup_ctl_endpoint) #ctlsocket to receive device successful startup signal

        self._ann_thread = threading.Thread(target=self._ann_loop) #start up announce thread
        self._ann_thread.start()

        
        #spawn task and result devices
        self._task_device = self._start_device('task', self.upstream_task_endpoint, self.downstream_task_endpoint, self._device_startup_ctl_endpoint, self.context) #task device
        self._result_device = self._start_device('result', self.upstream_result_endpoint, self.downstream_result_endpoint, self._device_startup_ctl_endpoint, self.context) #result device

        #get startup replies
        ann_ctlsocket.recv()
        device_ctlsocket.recv() #task device
        device_ctlsocket.recv() #result device

        ann_ctlsocket.close()
        device_ctlsocket.close()

        self._running = True
        self._shutdown_signalled = False

        log.debug('router: successful startup')

    def _ann_loop(self): 


        ctlsocket = self._make_signal_socket(self._announce_ctl_endpoint) #ctlsocket to receive shutdown signal

        # subscriber - a SUB socket to get announcements from server
        subscriber = self.context.socket(zmq.SUB)
        subscriber.setsockopt(zmq.SUBSCRIBE, '')
        subscriber.connect(self.upstream_announce_endpoint)

        # publisher - a PUB socket to send announcements to clients
        publisher = self.context.socket(zmq.PUB)
        publisher.bind(self.downstream_announce_endpoint)

        poller = zmq.Poller()
        poller.register(ctlsocket, zmq.POLLIN)
        poller.register(subscriber, zmq.POLLIN)

        self._signal_thread(self._startup_ctl_endpoint) # signal successful startup

        last_server_ping = None

        #'Work' loop
        try:
            while True:

                #TODO: add wait time - so that router shuts down if no response from server after specified time
                # for now, just simple receipt/publishing of announcements from server
                waittime = (self.server_heartbeat_interval*1000) #time to wait for server ping
                poll_results = dict(poller.poll(waittime))
                now = time.time()

                if poll_results.get(ctlsocket) == zmq.POLLIN:
                    if MSG_SHUTDOWN in recvall(ctlsocket):
                        publisher.send(MSG_SHUTDOWN)
                        log.debug('router: shutting down - forwarding message to clients')
                        return

                #Got announcement from server
                if poll_results.get(subscriber) == zmq.POLLIN:
                    announcements = recvall(subscriber)

                    if MSG_SHUTDOWN in announcements:
                        log.debug('router: received shutdown message from server')
                        self._shutdown() #Note - this function ends up sending shutdown message to clients
                        

                    if MSG_PING in announcements:
                        log.debug('router: pinged by server; forwarding ping to clients')
                        publisher.send(MSG_PING)
                        last_server_ping = now

                if last_server_ping is not None and (now - last_server_ping) >= 3*self.server_heartbeat_interval:
                    log.debug('Router: no ping from server; shutting down')
                    self._shutdown()

        #clean up
        finally:

            poller.unregister(ctlsocket)
            poller.unregister(subscriber)

            ctlsocket.close()
            subscriber.close()
            publisher.close()

            log.debug('router: exiting announce loop')


    def _shutdown(self): 

        if not self._shutdown_signalled:

            self._shutdown_signalled = True
            self._running = False

            

            #shut down devices
            self._task_device.shutdown()
            self._result_device.shutdown()

            #send shutdown message to announce loop
            self._signal_thread(self._announce_ctl_endpoint, MSG_SHUTDOWN)

    def _wait_for_shutdown(self):

        self._ann_thread.join()
        self._task_device.join()
        self._result_device.join()


    def _start_device(self, name, upstream_endpoint, downstream_endpoint, startup_ctl_endpoint, context):
        '''start up device'''

        device = ZMQDevice(name, upstream_endpoint, downstream_endpoint, startup_ctl_endpoint, context=context)
        device.startup()

        return device


    def shutdown(self):
        self.remove_router_info_file()
        self._shutdown()
        self._wait_for_shutdown()
        

    def __enter__(self):
        
        return self
        
    def __exit__(self, exc_type, exc_val, exc_traceback):
        self.remove_router_info_file()
        return False












        