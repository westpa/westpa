from __future__ import division, print_function; __metaclass__ = type

import sys, os, tempfile, logging, socket, multiprocessing, cPickle, threading, time, uuid, signal
import argparse
import collections
import zmq
from wemd.work_managers import WEMDWorkManager

from collections import deque

log = logging.getLogger(__name__)

DEFAULT_COMMAND_PORT = 23811
DEFAULT_STATUS_PORT  = 23812
DEFAULT_TASK_PORT    = 23813
DEFAULT_RESULTS_PORT = 23814

class Task:        
    def __init__(self, task_type, operand, task_id = None):
        self.task_id = task_id or uuid.uuid4()
        self.task_type = task_type
        self.operand = operand
                
    def __hash__(self):
        return hash(self.task_id)
        
    def __eq__(self, other):
        return self.task_id == other.task_id

    def __repr__(self):
        return '<Task 0x{id:x}: {self.task_id!s} {self.task_type!r} operand={self.operand!r}>'.format(id=id(self), self=self)
            
class ZMQBase:
    def __init__(self, zmq_context):
        self.context = zmq_context
        self.hostname = socket.gethostname()
        self.nodespec = '{}#{:d}:0x{:x}'.format(self.hostname, os.getpid(), threading.current_thread().ident)
        self.node_id = str(uuid.uuid4())
        
        # An inproc PUB socket for asynchronous inter-thread signaling WITHIN THIS CLASS
        self._intsig_endpoint = 'inproc://_intsig_{:x}'.format(id(self))
        self._intsig = self.context.socket(zmq.PUB)
        self._intsig.setsockopt(zmq.LINGER, 100)
        self._intsig.bind(self._intsig_endpoint)
                
    def listen_intsig(self, prefix=''):
        socket = self.context.socket(zmq.SUB)
        socket.setsockopt(zmq.SUBSCRIBE, prefix)
        socket.setsockopt(zmq.LINGER,0)
        socket.connect(self._intsig_endpoint)
        return socket
    
    def intsig(self, message):
        '''Send an in-process signal to all registered listeners'''
        self._intsig.send(message)
    
        
class ZWorkProvider(ZMQBase):
    def __init__(self, zmq_context, task_endpoint, results_endpoint, ann_tasks_endpoint):
        ZMQBase.__init__(self, zmq_context)
                
        self.task_endpoint = task_endpoint
        self.results_endpoint = results_endpoint
        self.ann_tasks_endpoint = ann_tasks_endpoint
        
        # create sockets in the threads in which they are used
        #self.task_socket = self.context.socket(zmq.REP)
        #self.results_socket = self.context.socket(zmq.PULL)
        #self.ann_tasks_socket = self.context.socket(zmq.PUB)
        #self.ann_tasks_socket.setsockopt(zmq.LINGER, 0)       
        #self.task_socket.bind(task_endpoint)
        #self.results_socket.bind(results_endpoint)
        #self.ann_tasks_socket.bind(ann_tasks_endpoint)
        
        self.task_socket = None
        self.results_socket = None
        self.ann_tasks_socket = None
        
        self.do_dispatch = True
        self.dispatch_thread = None
        
        self.do_collect = True
        self.collect_thread = None
                
        self.announce_interval = 100
        self.results_check_interval = 100

        # incoming tasks, undispatched
        self.task_queue = collections.deque()
        
        # Tasks (or rather, task_ids) which have been dispatched but not returned
        self.pending_tasks = set()
                
        # Mapping of task_id to results
        self.results = {}

    def handle_work_request(self):
        req = self.task_socket.recv()
        if not req.startswith('request'):
            return
        
        _reqtext, ntasktext = req.split()
        n_tasks = int(ntasktext)
        
        to_send = []
        while len(to_send) < n_tasks:
            try:
                # Atomic pops from task_queue
                to_send.append(self.task_queue.popleft())
            except IndexError:
                break        
        if to_send:
            log.debug('sending {:d} tasks: {!r}'.format(len(to_send), [task.task_id for task in to_send]))
            self.task_socket.send('tasks {} {}'.format(len(to_send), self.results_endpoint), zmq.SNDMORE)
            self.task_socket.send_pyobj(to_send)
            self.pending_tasks.update(task.task_id for task in to_send)
        else:
            log.debug('no tasks to send')
            self.task_socket.send('tasks 0', zmq.SNDMORE)
            self.task_socket.send('')
            
    def handle_results_receipt(self):
        results_queued = True
        while results_queued:
            try:
                task_id, results = self.results_socket.recv_pyobj(flags=zmq.NOBLOCK)
            except zmq.ZMQError as e:
                log.debug('no more results available: {}'.format(e))
                results_queued = False
            else:                
                log.debug('received results for task {!r}'.format(task_id))
                try:
                    self.pending_tasks.remove(task_id)
                except KeyError:
                    log.error('results received for unknown task {!r}; ignoring'.format(task_id))
                else:
                    # Atomic assignment
                    self.results[task_id] = results
                    # Signal to listeners to check for results
                    self.intsig('results {}'.format(task_id))        
        
    def wait(self, task_id):
        log.debug('waiting on {!r}'.format(task_id))
        try:
            results = self.results.pop(task_id)
            log.debug('wait complete for {!r}'.format(task_id))
            return results
        except KeyError:
            intsig = self.listen_intsig('results')
            try:
                while True:
                    try:
                        results = self.results.pop(task_id)
                        log.debug('wait complete for {!r}'.format(task_id))
                        return results
                    except KeyError:
                        # Not available - wait to hear about results                        
                        (rl,wl,el) = zmq.select([intsig], [], [], timeout=0.1)
                        if not rl:
                            log.debug('timeout waiting on {!r}'.format(task_id))
                            log.debug('results keys: {!r}'.format(list(self.results.keys())))
                        else:
                            assert rl == [intsig]
                            msg = intsig.recv()
            finally:
                intsig.close()
                pass
    
    def wait_all(self, task_ids):
        return [self.wait(task_id) for task_id in task_ids]
        
    def dispatch(self, task):    
        self.task_queue.append(task)
            
        if self.dispatch_thread is None:
            self.dispatch_thread = threading.Thread(target=self.dispatch_loop)
            self.dispatch_thread.daemon = True
            self.dispatch_thread.start()
            
    def dispatch_all(self, tasks):
        for task in tasks:
            # Atomic append
            self.task_queue.append(task)
            
        if self.dispatch_thread is None:
            self.dispatch_thread = threading.Thread(target=self.dispatch_loop)
            self.dispatch_thread.daemon = True
            self.dispatch_thread.start()            
            
    def collect(self):
        if self.collect_thread is None:
            self.collect_thread = threading.Thread(target=self.collect_loop)
            self.collect_thread.daemon = True
            self.collect_thread.start()
                    
    def shutdown(self):
        self.do_dispatch = False
        self.do_collect = False
        self.intsig('exit')
    
    def dispatch_loop(self):
        intsig = self.listen_intsig()
        self.ann_tasks_socket = self.context.socket(zmq.PUB)
        self.ann_tasks_socket.setsockopt(zmq.LINGER,0)
        self.ann_tasks_socket.bind(self.ann_tasks_endpoint)
        self.task_socket = self.context.socket(zmq.REP)
        self.task_socket.bind(self.task_endpoint)
            
        try:
            poller = zmq.Poller()
            poller.register(self.task_socket, zmq.POLLIN)
            poller.register(intsig, zmq.POLLIN)
            
            if self.task_queue:
                self.ann_tasks_socket.send('tasks available {}'.format(self.task_endpoint))
            
            while self.do_dispatch:
                poll_results = set(poller.poll(self.announce_interval))
                if (self.task_socket, zmq.POLLIN) in poll_results:
                    self.handle_work_request()
                if (intsig, zmq.POLLIN) in poll_results and intsig.recv() == 'exit':
                    return
                
                if not poll_results:
                    # timeout
                    log.debug('timeout')
                    if self.task_queue:
                        self.ann_tasks_socket.send('tasks available {}'.format(self.task_endpoint))
        finally:
            self.task_socket.close()
            self.ann_tasks_socket.close()
            del self.task_socket, self.ann_tasks_socket
            self.task_socket = None
            self.ann_tasks_socket = None
            self.dispatch_thread = None
            
    def collect_loop(self):
        intsig = self.listen_intsig()
        self.results_socket = self.context.socket(zmq.PULL)
        self.results_socket.bind(self.results_endpoint)

        try:        
            poller = zmq.Poller()
            poller.register(self.results_socket, zmq.POLLIN)
            poller.register(intsig, zmq.POLLIN)
            
            while self.do_collect:
                poll_results = set(poller.poll(1000))
                if (self.results_socket, zmq.POLLIN) in poll_results:
                    self.handle_results_receipt()
                    
                if (intsig, zmq.POLLIN) in poll_results and intsig.recv() == 'exit':
                    return
                
                if not poll_results:
                    log.debug('timeout')
        finally:
            self.results_socket.close()
            del self.results_socket
            self.results_socket = None
            self.collect_thread = None
            
class ZWorkConsumer(ZMQBase):
    def __init__(self, zmq_context, ann_in_endpoint, propagator):
        ZMQBase.__init__(self, zmq_context)
        
        self.propagator = propagator
        self.blocksize = 1
        
        self.ann_in_endpoint = ann_in_endpoint
        self.ann_in_socket = None
        
        #self.ann_in_socket = self.context.socket(zmq.SUB)
        #self.ann_in_socket.setsockopt(zmq.SUBSCRIBE,'')
        #self.ann_in_socket.connect(self.ann_in_endpoint)
                
        self.do_work = True
        self.work_thread = None
        
    def shutdown(self):
        self.do_work = False
        self.intsig('exit')
        
    def work(self):
        if self.work_thread is None:
            self.work_thread = threading.Thread(target=self.work_loop)
            self.work_thread.daemon = True
            self.work_thread.start()
    
    def retrieve_work(self, task_socket):
        task_socket.send('request {}'.format(self.blocksize))
        rpl = task_socket.recv_multipart()
        if not rpl[0].startswith('tasks'):
            log.error('invalid reply to task request; ignoring')
            return (None, None)
        
        taskfields = rpl[0].split()
        n_tasks = int(taskfields[1])
        if not n_tasks:
            log.debug('no tasks available')
            return (None,None)
        else:
            results_endpoint = taskfields[2]
        
        if n_tasks > self.blocksize:
            log.error('invalid number of tasks received; ignoring')
            return (None, None)
        
        payload = cPickle.loads(rpl[1])
        if len(payload) != n_tasks:
            log.error('payload contained wrong number of tasks; ignoring')
            return (None, None)
        
        return (results_endpoint, payload)
            
        
    def work_loop(self):
        self.ann_in_socket = self.context.socket(zmq.SUB)
        self.ann_in_socket.setsockopt(zmq.SUBSCRIBE,'')
        self.ann_in_socket.connect(self.ann_in_endpoint)
        intsig = self.listen_intsig()
        
        try:
        
            poller = zmq.Poller()
            poller.register(self.ann_in_socket, zmq.POLLIN)
            poller.register(intsig, zmq.POLLIN)
            
            while self.do_work:
                poll_results = set(poller.poll(1000))
                
                if not poll_results:
                    log.debug('timeout')
    
                
                if (intsig, zmq.POLLIN) in poll_results and intsig.recv() == 'exit':
                    self.work_thread = None
                    return
                
                if (self.ann_in_socket, zmq.POLLIN) in poll_results:
                    msg = self.ann_in_socket.recv()
                    if not msg.startswith('tasks available'):
                        log.error('invalid reply to task available query; ignoring') 
                        continue
                    _prefix, task_endpoint = msg.rsplit(None, 1)
                    
                    # Loop until we can get no more work from the given  provider, then sleep
                    # until the next announcement
                    task_socket = self.context.socket(zmq.REQ)
                    task_socket.connect(task_endpoint)
                    try:
                        results_endpoint, tasks = self.retrieve_work(task_socket)
                        while (results_endpoint, tasks) != (None, None):
                            # Connect the results socket (in case it fails)
                            #sys.stdout.write('{!r}'.format((results_endpoint,tasks)))
                            results_socket = self.context.socket(zmq.PUSH)
                            results_socket.connect(results_endpoint)
                            log.debug('received task(s) {!r}'.format([task.task_id for task in tasks]))
                            for task in tasks:
                                # Do work
                                log.debug('performing {!r} for {!r}'.format(task.task_type, task.task_id))
                                if task.task_type == 'propagate':
                                    segments = list(task.operand)
                                    self.propagator.propagate(segments)
                                    result = segments
                                    
                                # Ship out results one at a time - upstream expects a tuple of (task_id, result)
                                log.debug('returning results for {!r}'.format(task.task_id))
                                results_socket.send_pyobj((task.task_id, result))
                            
                            # And ask for more
                            results_endpoint, tasks = self.retrieve_work(task_socket)
                    finally:
                        task_socket.close()
        finally:
            self.ann_in_socket.close()
            del self.ann_in_socket
            self.ann_in_socket = None
            self.work_thread = None
                            
                
class ZMQWorkManager(WEMDWorkManager):
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        super(ZMQWorkManager,self).__init__(sim_manager)
        
        runtime_config = sim_manager.runtime_config
        self.mode = 'master'
        self.upstream_host  = runtime_config.get('zmq_work_manager.master', socket.gethostname())
        self.command_port = runtime_config.get_int('zmq_work_manager.command_port', DEFAULT_COMMAND_PORT)
        self.status_port  = runtime_config.get_int('zmq_work_manager.status_port', DEFAULT_STATUS_PORT)
        self.task_port    = runtime_config.get_int('zmq_work_manager.task_port', DEFAULT_TASK_PORT)
        self.results_port = runtime_config.get_int('zmq_work_manager.results_port', DEFAULT_RESULTS_PORT)
        self.n_workers = 1
        self.blocksize = 1
        
        self.context = None # zmq context
        self.command_endpoint = None # master commands pub/sub
        self.task_endpoint = None
        self.status_endpoint = None
        self.results_endpoint = None
        self.ann_results_endpoint = 'inproc://ann_results'
        self.ann_results_socket = None
        
        self.command_socket = None
        self.do_listen_commands = True
        self.command_thread = None
        
        self.zmaster = None
        self.zclient = None

    def parse_aux_args(self, aux_args, do_help = False, mode='master'):
        parser = argparse.ArgumentParser(usage='%(prog)s [NON_WORK_MANAGER_OPTIONS] [OPTIONS]',
                                         add_help=False)
        
        runtime_config = self.sim_manager.runtime_config
        parser.add_argument('-H', '--host', default=socket.gethostname(),
                            help='Upstream (master) host (default: %(default)s)')
        parser.add_argument('--cport', type=int, default=DEFAULT_COMMAND_PORT,
                            help='Port which the master uses to send commands to workers (default: %(default)s)')
        parser.add_argument('--sport', type=int, default=DEFAULT_STATUS_PORT,
                            help='Port on which status announcements are broadcast (default: %(default)s)')
        parser.add_argument('--tport', type=int, default=DEFAULT_TASK_PORT,
                            help='Port where workers contact the master to receive work (default: %(default)s)')
        parser.add_argument('--rport', type=int, default=DEFAULT_RESULTS_PORT,
                            help='Port where workers return completed tasks to the master (default: %(default)s)')
        parser.add_argument('-n', '-np', '-nw', type=int, dest='n_workers', default=multiprocessing.cpu_count(),
                            help='Number of worker processes to run on this host (default: %(default)s)')
        parser.add_argument('-b', '--blocksize', type=int, dest='blocksize', default=1,
                            help='Number of segments for each worker to run at once (default: %(default)s)')

        if do_help:
            parser.print_help()
            sys.exit(0)
        args, extra_args = parser.parse_known_args(aux_args)
        
        # zmq expects an IP address, not a hostname, for TCP transports
        # otherwise, the cryptic "no such device" error happens
        self.upstream_host =  runtime_config['zmq_work_manager.upstream_host'] = socket.gethostbyname(args.host)
        self.command_port =   runtime_config['zmq_work_manager.command_port'] = args.cport
        self.task_port =      runtime_config['zmq_work_manager.task_port'] = args.tport
        self.results_port =   runtime_config['zmq_work_manager.results_port'] = args.rport
        self.n_workers =      runtime_config['zmq_work_manager.n_workers'] = args.n_workers
        self.blocksize =      runtime_config['zmq_work_manager.blocksize'] = args.blocksize
        
        self.command_endpoint = 'tcp://{}:{}'.format(self.upstream_host, self.command_port)
        self.status_endpoint  = 'tcp://{}:{}'.format(self.upstream_host, self.status_port)
        self.task_endpoint    = 'tcp://{}:{}'.format(self.upstream_host, self.task_port)
        self.results_endpoint = 'tcp://{}:{}'.format(self.upstream_host, self.results_port)
        
        self.mode = mode
        
        return extra_args
        
    def shutdown(self, exit_code = 0):
        if self.mode == 'master':
            self.command_socket.send('cmd_shutdown')
            self.do_listen_commands = False
            self.zmaster.shutdown()
            
    def listen_commands_loop(self):
        '''Listen for and respond to commands from upstream'''
        #intsig = self.listen_intsig()
        poller = zmq.Poller()
        poller.register(self.command_socket, zmq.POLLIN)
        #poller.register(intsig, zmq.POLLIN)
        while self.do_listen_commands:
            poll_results = set(poller.poll(1000))
                        
            if (self.command_socket, zmq.POLLIN) in poll_results:
                command = self.command_socket.recv()
                if command == 'cmd_shutdown':
                    self.zworker.shutdown()
                    return
                
            if not poll_results:
                log.debug('timeout')

            
    def listen_commands(self):
        if self.command_thread is None:
            self.command_thread = threading.Thread(target=self.listen_commands_loop)
            self.command_thread.daemon = True
            self.command_thread.start()
            
    def run_worker(self):
        assert self.zworker is not None
        assert self.mode == 'worker'
        self.listen_commands()
        self.zworker.work()
        return self.zworker.work_thread.join()
    
    def prepare(self):
        if self.mode == 'master':
            # Run workers on this machine, and also run as master
            for iworker in xrange(0, self.n_workers):
                log.debug('spawning worker {} of {}'.format(iworker+1, self.n_workers))
                pid = os.fork()
                if pid == 0:
                    # in child
                    self.mode = 'worker'
                    self.context = zmq.Context()
                    self.command_socket = self.context.socket(zmq.SUB)
                    self.command_socket.setsockopt(zmq.SUBSCRIBE, '')
                    self.command_socket.connect(self.command_endpoint)
                    self.zworker = ZWorkConsumer(self.context, self.status_endpoint, self.sim_manager.propagator)
                    log.debug('worker {} is PID {}'.format(iworker+1, os.getpid()))
                    return
                            
            # Run master
            assert self.context is None
            log.debug('master is PID {}'.format(os.getpid()))
            self.context = zmq.Context()
            self.command_socket = self.context.socket(zmq.PUB)
            self.command_socket.setsockopt(zmq.LINGER, 0)
            self.command_socket.bind(self.command_endpoint)
            self.ann_results_socket = self.context.socket(zmq.SUB)
            self.ann_results_socket.bind(self.ann_results_endpoint)
            self.zmaster = ZWorkProvider(self.context, self.task_endpoint, self.results_endpoint, self.status_endpoint)
            self.zmaster.collect()
        else:
            raise NotImplementedError
    
    def propagate(self, segments):
        assert self.zmaster.collect_thread is not None
        blocks = [segments[i:i+self.blocksize] for i in xrange(0,len(segments),self.blocksize)]
        log.debug('dispatching {} segment(s) in {} block(s)'.format(len(segments), len(blocks)))
        tasks = [Task('propagate', block) # task_id=('propagate', tuple((segment.n_iter, segment.seg_id) for segment in block))) 
                 for block in blocks]
        self.zmaster.dispatch_all(tasks)
        results = self.zmaster.wait_all(task.task_id for task in tasks)
        incoming_segments = []
        for block in results:
            incoming_segments.extend(block)
        
        outgoing_by_id = {segment.seg_id: segment for segment in segments}
        incoming_by_id = {segment.seg_id: segment for segment in incoming_segments}
        assert set(outgoing_by_id) == set(incoming_by_id)
        
        for incoming_segment in incoming_segments:
            orig_segment = outgoing_by_id[incoming_segment.seg_id]
            orig_segment.status = incoming_segment.status
            orig_segment.walltime = incoming_segment.walltime
            orig_segment.cputime  = incoming_segment.cputime
            orig_segment.pcoord[...] = incoming_segment.pcoord[...]
