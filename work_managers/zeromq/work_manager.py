'''
The ZeroMQ work manager, which combines server and client functionality.
'''

from __future__ import division, print_function; __metaclass__ = type


import os, sys, socket, uuid, multiprocessing, json, re, atexit, logging, warnings

import work_managers
from work_managers import WorkManager

from server import ZMQServer
from client import ZMQClient

from . import (DEFAULT_HANGCHECK_INTERVAL, DEFAULT_MAX_TASKQUEUE_SIZE, DEFAULT_SERVER_HEARTBEAT_INTERVAL, 
               DEFAULT_SHUTDOWN_TIMEOUT, DEFAULT_TASK_TIMEOUT, DEFAULT_TASKQUEUE_WAIT)

log = logging.getLogger(__name__)

class ZMQWorkManager(ZMQServer,WorkManager):
    write_server_info = True
    
    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 

        wm_group = parser.add_argument_group('options for ZeroMQ ("zmq") work manager or router')
        wm_group.add_argument(wmenv.arg_flag('zmq_mode'), metavar='MODE', choices=('server', 'client'),
                              help='Operate as a server (MODE=server) or a client (MODE=client).')
        wm_group.add_argument(wmenv.arg_flag('zmq_info'), metavar='INFO_FILE',
                              help='Store server information in INFO_FILE. (specify for server or routers only)'
                                   'This is helpful if running server and clients or routers on multiple '
                                   'machines which share a filesystem, as explicit hostnames/ports are not required')
        wm_group.add_argument(wmenv.arg_flag('zmq_task_endpoint'), metavar='TASK_ENDPOINT',
                              help='''Bind server to given ZeroMQ endpoint to distribute tasks downstream (specify for servers or routers only).''')
        wm_group.add_argument(wmenv.arg_flag('zmq_result_endpoint'), metavar='RESULT_ENDPOINT',
                              help='''Bind server to given ZeroMQ endpoint to receive results from downstream (specify for servers or routers only).''')
        wm_group.add_argument(wmenv.arg_flag('zmq_announce_endpoint'), metavar='ANNOUNCE_ENDPOINT',
                              help='''Bind server to given ZeroMQ endpoint to send anouncements downstream (specify for servers or routers only).''')
        wm_group.add_argument(wmenv.arg_flag('zmq_heartbeat_interval'), metavar='INTERVAL',
                              help='''If a client or router has not
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
         
        if wmenv.get_val('zmq_mode','server').lower() == 'client':
            return ZMQClient.from_environ()
        
        n_workers = wmenv.get_val('n_workers', multiprocessing.cpu_count(), int)
        task_timeout = wmenv.get_val('zmq_task_timeout', DEFAULT_TASK_TIMEOUT, int)
        heartbeat_interval = wmenv.get_val('zmq_heartbeat_interval', DEFAULT_SERVER_HEARTBEAT_INTERVAL, int)
        shutdown_timeout = wmenv.get_val('zmq_worker_shutdown_timeout', DEFAULT_SHUTDOWN_TIMEOUT, int)
        hangcheck_interval = wmenv.get_val('zmq_hangcheck_interval', DEFAULT_HANGCHECK_INTERVAL, int)
        client_comm_mode = wmenv.get_val('zmq_client_comm_mode')
        server_info_filename = wmenv.get_val('zmq_info') or wmenv.get_val('zmq_write_info', 'zmq_server_info_{}.json'.format(uuid.uuid4().hex))

        # if individual endpoints are named, we use these
        tests_old = [not bool(wmenv.get_val('zmq_task_endpoint')),
                 not bool(wmenv.get_val('zmq_result_endpoint')),
                 not bool(wmenv.get_val('zmq_announce_endpoint'))]
        tests_new = [not bool(wmenv.get_val('zmq_downstream_task_endpoint')),
                 not bool(wmenv.get_val('zmq_downstream_result_endpoint')),
                 not bool(wmenv.get_val('zmq_downstream_announce_endpoint'))]

        if all(tests_old) and all(tests_new):
            # Choose random ports
            task_endpoint = cls.canonicalize_endpoint('tcp://*')
            result_endpoint = cls.canonicalize_endpoint('tcp://*')
            announce_endpoint = cls.canonicalize_endpoint('tcp://*')
        elif (not all(tests_old) and any(tests_old)) or (not all(tests_new) and any(tests_new)):
            raise ValueError('either none or all three server endpoints must be specified')
        #Use new-style, unambiguous endpoint args
        elif not all(tests_new):
            task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_task_endpoint'))
            result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_result_endpoint'))
            announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_downstream_announce_endpoint'))
        else:
            task_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_task_endpoint'))
            result_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_result_endpoint'))
            announce_endpoint = cls.canonicalize_endpoint(wmenv.get_val('zmq_announce_endpoint'))
            
        return cls(n_workers, task_endpoint, result_endpoint, announce_endpoint, server_info_filename = server_info_filename,
                   server_heartbeat_interval=heartbeat_interval, task_timeout = task_timeout,
                   client_comm_mode = client_comm_mode, client_hangcheck_interval = hangcheck_interval,
                   worker_shutdown_timeout = shutdown_timeout)
        
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
                 write_server_info = True, server_info_filename=None, server_heartbeat_interval=None, 
                 task_timeout = None,
                 client_comm_mode = None, worker_shutdown_timeout=None, client_hangcheck_interval=None,
                 ):

        WorkManager.__init__(self)

        # 0 is permissible here, meaning dedicated server
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()
        self.n_workers = n_workers
        
        server_heartbeat_interval = server_heartbeat_interval or DEFAULT_SERVER_HEARTBEAT_INTERVAL
        task_timeout = task_timeout or DEFAULT_TASK_TIMEOUT
        worker_shutdown_timeout = worker_shutdown_timeout or DEFAULT_SHUTDOWN_TIMEOUT
        client_comm_mode = client_comm_mode or 'ipc'
        if not client_hangcheck_interval:
            client_hangcheck_interval = min(server_heartbeat_interval, task_timeout or sys.maxint)
        
        
        self.internal_client = None
        
        argtests = [master_task_endpoint is None, master_result_endpoint is None, master_announce_endpoint is None]
        if any(argtests) and not all(argtests):
            raise ValueError('endpoints must all be either specified or None (not mixed)')
        else:
            assign_endpoints = all(argtests)
            
        if assign_endpoints:
            master_task_endpoint = self.make_ipc_endpoint()
            master_result_endpoint = self.make_ipc_endpoint()
            master_announce_endpoint = self.make_ipc_endpoint()

        ZMQServer.__init__(self, master_task_endpoint, master_result_endpoint, master_announce_endpoint,
                             server_heartbeat_interval)            
        
        if n_workers > 0:
            # this node is both a master and a client; start workers
            self.internal_client = ZMQClient(master_task_endpoint, master_result_endpoint, master_announce_endpoint, 
                                             server_heartbeat_interval = server_heartbeat_interval,
                                             comm_mode = client_comm_mode,
                                             n_workers = n_workers, 
                                             task_timeout = task_timeout,
                                             shutdown_timeout = worker_shutdown_timeout,
                                             hangcheck_interval = client_hangcheck_interval)            
        
        if write_server_info:
            self.server_info_filename = server_info_filename or 'zmq_server_info_{}.json'.format(uuid.uuid4().hex)
            self.write_server_info(self.server_info_filename)
            atexit.register(self.remove_server_info_file)
        else:
            self.server_info_filename = None
                
                    
    def startup(self):
        if not self.running:
            self.running = True
            ZMQServer.startup(self)
            if self.internal_client is not None:
                self.internal_client.startup()
            
    def shutdown(self):
        if self.running:
            if self.internal_client is not None:
                self.internal_client.shutdown()
            ZMQServer.shutdown(self)
            self.remove_server_info_file()
            self.running = False
