from __future__ import division, print_function; __metaclass__ = type

import sys, os, logging, socket, multiprocessing, threading, time, traceback, signal, random
import argparse
from collections import deque
import zmq
import wemd
from wemd.work_managers import WEMDWorkManager, WMFuture

log = logging.getLogger(__name__)

def recvall(socket):
    messages = []
    while True:
        try:
            messages.append(socket.recv(flags=zmq.NOBLOCK))
        except zmq.ZMQError as err:
            if err.errno == zmq.EAGAIN:
                return messages
            else:
                raise
        
def recvall_pyobjs(socket):
    objs = []
    while True:
        try:
            objs.append(socket.recv_pyobj(flags=zmq.NOBLOCK))
        except zmq.ZMQError as err:
            if err.errno == zmq.EAGAIN:
                return objs
            else:
                raise

DEFAULT_ANN_PORT      = 23811 # Upstream-to-downstream announcements
DEFAULT_TASK_PORT     = 23812 # Task and result distribution
    
class ZMQWorkManager(WEMDWorkManager):
            
    def __init__(self):
        super(ZMQWorkManager,self).__init__()
         
        # Upstream hostname and IP
        self.upstream_hostname  = None
        self.upstream_host      = None
        
        # Ports
        self.ann_port           = None
        self.task_port          = None
        
        # (upstream) endpoints
        self.ann_endpoint       = None
        self.task_endpoint      = None
                        
        # ZeroMQ context
        self.context            = None
        
        # Information about this host
        self.this_hostname      = socket.gethostname()
        
        # Spawned processes
        self.spawned_pids       = list()
        
        # Thread-local storage
        self._tls = threading.local()
        
    def get_id_str(self):
        return '{:s}-{:d}'.format(self.this_hostname, os.getpid())
    
    id_str = property(get_id_str,None,None,None)
        
    def bind_thread_ctl(self, endpoint):
        log.debug('binding new socket for {}'.format(endpoint))
        ctlsocket = self.context.socket(zmq.PULL)
        ctlsocket.bind(endpoint)
        return ctlsocket
    
    def signal_thread(self, endpoint, message=''):
        try:
            ctlsocket = self._tls.signal_sockets[endpoint]
        except AttributeError:
            #log.debug('connecting new socket for {}'.format(endpoint))
            ctlsocket = self.context.socket(zmq.PUSH)
            #ctlsocket.setsockopt(zmq.HWM,1)
            ctlsocket.connect(endpoint)
            self._tls.signal_sockets = {endpoint: ctlsocket}
        except KeyError:
            #log.debug('connecting new socket for {}'.format(endpoint))            
            ctlsocket = self.context.socket(zmq.PUSH)
            #ctlsocket.setsockopt(zmq.HWM,1)
            ctlsocket.connect(endpoint)
            self._tls.signal_sockets[endpoint] = ctlsocket
        else:
            #log.debug('using already-open socket')
            pass
        
        ctlsocket.send(message)
                        
    def parse_aux_args(self, aux_args, do_help = False):
        parser = argparse.ArgumentParser(usage='%(prog)s [NON_WORK_MANAGER_OPTIONS] [OPTIONS]',
                                         add_help=False)
        
        group = parser.add_argument_group('zmq work manager options')

        mgroup = group.add_mutually_exclusive_group()
        mgroup.add_argument('--master', dest='mode', action='store_const', const=self.MODE_MASTER,
                            help='Run as a master process (which generates and distributes computational tasks; default).')
        mgroup.add_argument('--worker', dest='mode', action='store_const', const=self.MODE_WORKER,
                            help='Run as a worker process')
                        
        group.add_argument('-n',  type=int, dest='n_workers', default=multiprocessing.cpu_count(),
                           help='Number of worker processes to run on this host. Use 0 for a dedicated server '
                               +' or forwarder process. (Default: %(default)s)')
        group.add_argument('-H', '--host', default=wemd.rc.config.get('work_manager.zmq.master', socket.gethostname()),
                           help='Upstream (master/coordinator) host (default: %(default)s)')
        group.add_argument('--annport', type=int, 
                           default=wemd.rc.config.get_int('work_manager.zmq.ann_port', DEFAULT_ANN_PORT),
                           help='Port on which master makes announcements (default: %(default)s)')
        group.add_argument('--taskport', type=int, 
                           default=wemd.rc.config.get_int('work_manager.zmq.task_port', DEFAULT_TASK_PORT),
                           help='Port on master which remote workers contact to receive work (default: %(default)s)')
        
        parser.set_defaults(mode=self.MODE_MASTER)

        if do_help:
            parser.print_help()
            sys.exit(0)
        args, extra_args = parser.parse_known_args(aux_args)
        
        self.mode              = args.mode
        self.upstream_hostname = args.host
        self.upstream_host     = socket.gethostbyname(args.host)
        self.ann_port          = args.annport
        self.task_port         = args.taskport
        self.n_workers         = args.n_workers
        
        self.ann_endpoint = 'tcp://{:s}:{:d}'.format(self.upstream_host, self.ann_port)
        self.task_endpoint = 'tcp://{:s}:{:d}'.format(self.upstream_host, self.task_port)
        
        return extra_args
    
    def startup(self):
        # do all the forking here, and reassign either the master or worker to self.__class__ prior to
        # invoking start_threads()
        
        if self.mode == self.MODE_MASTER:
            if self.n_workers:
                # Not a dedicated master, so start workers
                self.start_workers()
                if self.mode == self.MODE_WORKER: 
                    return                
            self.make_master()
            self.start_threads()
            log.info('pid {:d} is master'.format(os.getpid()))
        elif self.mode == self.MODE_WORKER:
            self.start_workers()
            self.wait_workers()
        else:
            raise AssertionError('invalid mode {!r}'.format(self.mode))
        
        log.debug('pid {:d} has mode {!r}'.format(os.getpid(), self.mode))
        return self.mode
                
    def start_workers(self):
        for n in xrange(self.n_workers):
            log.info('forking worker process {:d} of {:d}'.format(n+1,self.n_workers))
            pid = os.fork()
            if not pid:
                log.info('pid {:d} is worker'.format(os.getpid()))
                
                # (re)initialize random number generator in this process
                random.seed()
                
                self.make_worker()
                self.start_threads()
                self.spawned_pids = []
                break
            else:
                self.spawned_pids.append(pid)
    
    def wait_workers(self):
        for pid in self.spawned_pids:
            os.waitpid(pid, 0)
    
    # These functions are defined in subclasses, which we magically transmute to
    # at the appropriate time
    def shutdown(self, exit_code=0):
        raise NotImplementedError
    
    def start_threads(self):
        raise NotImplementedError
    
    # Transmutation -- think of these as subsidiary __init__ functions
    def make_master(self):
        self.__class__ = ZMQMasterWorkManager
        self.task_queue = deque()
        self.pending_futures = dict()
        self.dr_ctl_endpoint = 'inproc://dr_ctl_{:x}'.format(id(self))
        self.ann_ctl_endpoint = 'inproc://ann_ctl_{:x}'.format(id(self))
        self.dr_thread = None
        self.ann_thread = None
        self.client_avail = False
        self.client_avail_check = 0.01
        self.work_avail_interval = 2
        self.context = zmq.Context.instance()
        self.mode = self.MODE_MASTER
        self.shutdown_called = False
        
    def make_worker(self):
        self.__class__ = ZMQWorkerWorkManager
        sys.stdin.close() # so we don't get Ctrl-C
        signal.signal(signal.SIGINT, signal.SIG_IGN) # so we can send SIGINT to a process group we lead
        os.setpgid(0,0) # make all children exist within the same process group
        self.rdp_ctl_endpoint = 'inproc://rdp_ctl_{:x}'.format(id(self))
        self.rdp_thread = None
        self.al_thread = None
        self.context = zmq.Context.instance()
        self.mode = self.MODE_WORKER
        self.shutdown_called = False
                
class ZMQMasterWorkManager(ZMQWorkManager):
    
    def start_threads(self):
        self.ann_thread = threading.Thread(target=self.announcer,name='zmq_announce')
        self.dr_thread = threading.Thread(target=self.dispatch_receive,name='zmq_dr')
        self.ann_thread.start()
        time.sleep(0.01)
        self.dr_thread.start()
                    
    def announcer(self):
        ann_ctl = self.bind_thread_ctl(self.ann_ctl_endpoint)
        ann_socket = self.context.socket(zmq.PUB)
        ann_socket.bind(self.ann_endpoint)
        poller = zmq.Poller()
        poller.register(ann_ctl,zmq.POLLIN)
        
        timeout = self.client_avail_check
        while not self.shutdown_called:
            poll_results = {fd: _flag for (fd, _flag) in poller.poll(timeout*1000)}
            
            if poll_results.get(ann_ctl) == zmq.POLLIN:
                messages = recvall(ann_ctl)
                if 'shutdown' in messages:
                    log.debug('shutdown received for master')
                    self.shutdown_called = True
                    if self.client_avail:
                        log.debug('sending shutdown signal to clients')
                        ann_socket.send('shutdown')
                    break
                    
            # The only other use for this is to announce work available, so just
            # wake up and send work available if there is any work and there are
            # any clients
            if self.client_avail:
                timeout = self.work_avail_interval
                if self.task_queue:
                    log.debug('sending tasks available signal')
                    ann_socket.send('task_avail')
        log.debug('exiting announcer')
            

    def dispatch_receive(self):        
        dr_ctl = self.bind_thread_ctl(self.dr_ctl_endpoint)
        
        task_socket = self.context.socket(zmq.REP)
        task_socket.bind(self.task_endpoint)
        
        poller = zmq.Poller()
        poller.register(dr_ctl, zmq.POLLIN)
        poller.register(task_socket, zmq.POLLIN)
        while not self.shutdown_called:
            poll_results = {fd: _flag for (fd, _flag) in poller.poll()}

            if poll_results.get(dr_ctl) == zmq.POLLIN:
                messages = recvall(dr_ctl)
                for message in messages:
                    if message == 'shutdown':
                        break
                del messages
                    
            if poll_results.get(task_socket) == zmq.POLLIN and not self.shutdown_called:
                ts_tuple = task_socket.recv_pyobj()
                
                message = ts_tuple[0]
                sender  = ts_tuple[1]
                payload = ts_tuple[2]
                
                if message == 'task':
                    log.debug('task request received from {!r}'.format(sender))
                    try:
                        task = self.task_queue.popleft()
                        task_id = task[0]
                    except IndexError:
                        log.debug('no task to send to {!r}'.format(sender))
                        task_socket.send_pyobj(None)
                    else:
                        log.debug('dispatching {!r} to {!r}'.format(task_id, sender))
                        task_socket.send_pyobj(task)
                        self.work_last_dispatched = time.time()
                        del task, task_id
                elif message == 'result':
                    task_socket.send('ack')
                    task_id, result_type = payload[0:2]
                    log.debug('result/exception received for task {!r} from {!r}'.format(task_id, sender))
                    ft = self.pending_futures.pop(task_id)
                    if result_type == 'result':
                        ft._set_result(payload[2])
                    elif result_type == 'exception':
                        ft._set_exception(*payload[2])
                    del task_id, result_type, ft
                elif message == 'ping':
                    log.debug('received ping from {!r}'.format(sender))
                    task_socket.send('ack')
                    self.client_avail = True
                    self.signal_thread(self.ann_ctl_endpoint)
                else:
                    log.error('unknown message received from {!s}: {!r}'.format(sender, message))
                    
                del ts_tuple, message, sender, payload
        task_socket.close(0)
        log.debug('exiting dispatch_receive')
                
    def submit(self, fn, *args, **kwargs):
        ft = WMFuture()
        task_id = ft.task_id
        tasktuple = (task_id, fn, args, kwargs)
        self.pending_futures[task_id] = ft
        self.task_queue.append(tasktuple)
        self.signal_thread(self.ann_ctl_endpoint)
        return ft
    
    def shutdown(self, exit_code=0):
        log.debug('signalling shutdown')
        self.signal_thread(self.dr_ctl_endpoint, 'shutdown')
        self.signal_thread(self.ann_ctl_endpoint, 'shutdown')
        self.wait_workers()

class ZMQWorkerWorkManager(ZMQWorkManager):
    
    def start_threads(self):
        self.rdp_thread = threading.Thread(target=self.request_perform_dispatch, name='zmq_rdp')
        self.la_thread = threading.Thread(target=self.listen_announcements, name='zmq_la')
        self.rdp_thread.start()
        self.la_thread.start()
        self.la_thread.join()
        self.rdp_thread.join()
        
    def interrupt_children(self):
        pgid = os.getpgid(0)
        log.debug('sending SIGINT to process group {:d}'.format(pgid))
        # Kill child processes
        os.killpg(pgid, signal.SIGINT)
        
    def listen_announcements(self):
        ann_socket = self.context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')
        ann_socket.connect(self.ann_endpoint)
        poller = zmq.Poller()
        poller.register(ann_socket,zmq.POLLIN)

        # Before doing anything else, contact the master and announce our presence
        # This will block until the master comes up
        task_socket = self.context.socket(zmq.REQ)
        task_socket.connect(self.task_endpoint)
        task_socket.send_pyobj(('ping', self.id_str, None))
        task_socket.recv()
        task_socket.close()
        del task_socket
        
        # Now that we have opened our announcement mailbox and informed the master we are here, we can start to
        # listen for work                
        while not self.shutdown_called:
            poll_results = {fd: _flag for (fd, _flag) in poller.poll()}
            
            if poll_results.get(ann_socket) == zmq.POLLIN:
                announcements = set(recvall(ann_socket))
                log.debug('received {:d} announcements'.format(len(announcements)))
                if 'shutdown' in announcements:
                    log.debug('shutdown message received at client')
                    self.interrupt_children()
                    self.shutdown_called = True
                    self.signal_thread(self.rdp_ctl_endpoint, 'shutdown')
                    break
                elif 'task_avail' in announcements:
                    self.signal_thread(self.rdp_ctl_endpoint, 'task_avail')
        log.debug('exiting listen_announcements')
    
    def request_perform_dispatch(self):
        rdp_ctl_socket = self.bind_thread_ctl(self.rdp_ctl_endpoint)
        
        poller = zmq.Poller()
        poller.register(rdp_ctl_socket,zmq.POLLIN)
                
        while not self.shutdown_called:
            poll_results = {fd: _flag for (fd, _flag) in poller.poll()}
                          
            if poll_results.get(rdp_ctl_socket) == zmq.POLLIN:
                messages = recvall(rdp_ctl_socket)
                if 'shutdown' in messages:
                    log.debug('shutdown message received in request_perform_dispatch')
                    return
                elif 'task_avail' in messages:
                    log.debug('task available message received')
                    task_avail = True            
                    while task_avail and not self.shutdown_called:
                        task_socket = self.context.socket(zmq.REQ)
                        task_socket.connect(self.task_endpoint)
                        try:
                            log.debug('sending request for task')                        
                            task_socket.send_pyobj(('task', self.id_str, None))
                            task = task_socket.recv_pyobj()
                            task_socket.close()
                        finally:
                            del task_socket
                            
                        if task is None:
                            log.debug('no further task available')
                            task_avail = False
                        else:
                            log.debug('task received: {!r}'.format(task))
                            (task_id, fn, args, kwargs) = task
                            try:
                                result = fn(*args, **kwargs)
                            except Exception as exception:
                                result_tuple = (task_id, 'exception', (exception, traceback.format_exc()))
                            else:
                                result_tuple = (task_id, 'result', result)
                            
                            # if we have been told to shut down in the meantime, don't try to send a 
                            # result (that will hang)
                            if not self.shutdown_called:
                                result_socket = self.context.socket(zmq.REQ)
                                result_socket.connect(self.task_endpoint)
                                try:
                                    log.debug('dispatching result')
                                    result_socket.send_pyobj(('result', self.id_str, result_tuple))
                                    result_socket.recv()
                                    result_socket.close()
                                finally:
                                    del result_socket
                del messages
        log.debug('exiting request_perform_dispatch')
                            
    def shutdown(self, exit_code):
        return 
