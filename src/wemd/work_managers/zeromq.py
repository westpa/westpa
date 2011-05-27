from __future__ import division, print_function; __metaclass__ = type

import sys, os, tempfile, logging, socket, multiprocessing, cPickle, threading, time, uuid, signal
import argparse
import collections
import zmq
from wemd.work_managers import WEMDWorkManager

log = logging.getLogger(__name__)

DEFAULT_ANN_PORT     = 23811
DEFAULT_TASK_PORT    = 23812
DEFAULT_RESULTS_PORT = 23813

def assert_flushed(socket):
    try:
        socket.recv(flags=zmq.NOBLOCK)
    except zmq.ZMQError as zme:
        if zme.errno != zmq.EAGAIN:
            raise
    else:
        raise AssertionError('socket has data remaining')

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
    def __init__(self, context=None):
        self.context = context
        self._ipc_endpoints = []

    def is_valid_ipc_endpoint(self, endpoint):
        if not endpoint.startswith('ipc://'):
            return False
        else:
            return os.path.exists(endpoint[6:])

    def make_ipc_endpoint(self):
        (fd, socket_path) = tempfile.mkstemp()
        os.close(fd)
        endpoint = 'ipc://{}'.format(socket_path)
        self._ipc_endpoints.append(endpoint)
        return endpoint
        
    def remove_ipc_endpoints(self):
        for endpoint in self._ipc_endpoints:
            assert endpoint.startswith('ipc://')
            socket_path = endpoint[6:]  
            try:
                os.unlink(socket_path)
            except OSError as e:
                log.debug('could not unlink IPC endpoint {!r}: {}'.format(socket_path, e))
            else:
                log.debug('unlinked IPC endpoint {!r}'.format(socket_path))
                
    def flush_socket(self, socket):
        messages = []
        while True:
            try:
                messages.append(socket.recv(flags=zmq.NOBLOCK))
            except zmq.ZMQError as zme:
                if zme.errno == zmq.EAGAIN:
                    return messages
                else:
                    raise

class Master(ZMQBase):
    def __init__(self, zmq_context, task_endpoint, results_endpoint, ann_endpoint):
        ZMQBase.__init__(self, zmq_context)
                
        self.task_endpoint = task_endpoint
        self.results_endpoint = results_endpoint
        self.ann_endpoint = ann_endpoint
        
        # How often do we broadcast that we have work available, to handle slow starting clients (in s)
        self.announce_interval = 5.0
                                
        # How often do we check for results (in s)
        self.results_check_interval = 0.05
        
        self.do_main_loop = True
        self.main_loop_thread = None
        self._work_avail_endpoint = 'inproc://_work_avail_{:x}'.format(id(self))
        self._results_avail_endpoint = 'inproc://_results_avail_{:x}'.format(id(self))
        self._main_wakeup_endpoint = 'inproc://_wakeup_main_{:x}'.format(id(self))

        # incoming tasks, undispatched
        self.task_queue = collections.deque()
                        
        # Mapping of task_id to results (i.e. returned results)
        # the main loop inserts items, and wait() removes them
        self.results = {}
        
    def main_loop(self):
        log.debug('starting main loop')
        ann_socket = self.context.socket(zmq.PUB)
        results_socket = self.context.socket(zmq.PULL)
        results_avail_socket = self.context.socket(zmq.PUB)
        work_avail_socket = self.context.socket(zmq.SUB)
        work_avail_socket.setsockopt(zmq.SUBSCRIBE,'')
        task_socket = self.context.socket(zmq.REP)
        wakeup_socket = self.context.socket(zmq.SUB)
        wakeup_socket.setsockopt(zmq.SUBSCRIBE, '')
        
        ann_socket.bind(self.ann_endpoint)
        results_socket.bind(self.results_endpoint)
        task_socket.bind(self.task_endpoint)
        wakeup_socket.bind(self._main_wakeup_endpoint)
        results_avail_socket.bind(self._results_avail_endpoint)
        work_avail_socket.bind(self._work_avail_endpoint)
                
        poller = zmq.Poller()
        poller.register(results_socket, zmq.POLLIN)
        poller.register(task_socket, zmq.POLLIN)
        poller.register(wakeup_socket, zmq.POLLIN)
        poller.register(work_avail_socket, zmq.POLLIN)
                
        log.debug('waiting for sockets to settle')
        time.sleep(1)        
        
        last_announcement = 0
        while self.do_main_loop:            
            # Wait on messages from clients (or the shutdown signal)
            poll_results = set(fd for (fd,_flag) in poller.poll(self.announce_interval*1000))
            
            if wakeup_socket in poll_results:
                # Just flush, and make sure we re-evaluate self.do_main_loop and announce work if it's time
                log.debug('wakeup!')
                wakeup_socket.recv()

            if results_socket in poll_results:
                log.debug('receiving results')
                task_id, results = results_socket.recv_pyobj()
                log.debug('received results for task {!r}'.format(task_id))
                #log.debug('  results: {!r}'.format(results))
                
                # Atomic assignment
                self.results[task_id] = results
                # Notify watchers that new results are available
                results_avail_socket.send('')

            if task_socket in poll_results:
                log.debug('receiving task request')
                req = task_socket.recv()
                assert req.startswith('request')
                
                _reqtext, ntasktext = req.split()
                n_tasks = int(ntasktext)
                
                to_send = []
                while len(to_send) < n_tasks:
                    try:
                        # Atomic pops from task_queue
                        to_send.append(self.task_queue.popleft())
                    except IndexError:
                        break
                log.debug('sending {} task(s)'.format(len(to_send)))
                task_socket.send_pyobj((self.results_endpoint, to_send))
            
            now = time.time()
            if self.task_queue and (work_avail_socket in poll_results or now - last_announcement >= self.announce_interval):
                if work_avail_socket in poll_results:
                    work_avail_socket.recv()
                log.debug('announcing tasks available')
                ann_socket.send('tasks available {}'.format(self.task_endpoint))
                last_announcement = now
            
        else:
            log.debug('sending termination signal')
            ann_socket.send('shutdown')

            
    def wait(self, task_id):        
        # First, see if the result is already waiting for us
        try:
            result = self.results[task_id]
        except KeyError:
            pass
        else:
            del self.results[task_id]
            log.debug('result ready immediately for task {!r}'.format(task_id))
            return result
            
        # If the result is not waiting, we wait on it
        results_avail = self.context.socket(zmq.SUB)
        results_avail.setsockopt(zmq.SUBSCRIBE, '')
        results_avail.connect(self._results_avail_endpoint)
        
        poller = zmq.Poller()
        poller.register(results_avail, zmq.POLLIN)
        while True:
            poll_results = set(fd for (fd, _flags) in poller.poll(self.results_check_interval*1000))
            
            if results_avail in poll_results:
                results_avail.recv()
                #self.flush_socket(results_avail)
                            
            try:
                result = self.results[task_id]
            except KeyError:
                log.debug('result not ready for task {!r}'.format(task_id))
                continue
            else:
                log.debug('result ready for task {!r}'.format(task_id))
                del self.results[task_id]
                return result
        
    def wait_all(self, task_ids):
        return [self.wait(task_id) for task_id in list(task_ids)]
        
    def start_main_thread(self):
        if self.main_loop_thread is None:
            self.do_main_loop = True
            self.main_loop_thread = threading.Thread(target=self.main_loop)
            self.main_loop_thread.daemon = False
            self.main_loop_thread.start()
                
    def dispatch(self, task=None):
        self.start_main_thread()        
        if task is not None:
            self.task_queue.append(task)
            socket = self.context.socket(zmq.PUB)
            socket.connect(self._work_avail_endpoint)
            socket.send('')                    
            
    def dispatch_all(self, tasks):
        self.start_main_thread()
        for task in tasks:
            # Atomic append; extend() is probably not a good idea here
            self.task_queue.append(task)
        
        if tasks:
            socket = self.context.socket(zmq.PUB)
            socket.connect(self._work_avail_endpoint)
            socket.send('')
                            
    def shutdown(self):
        log.debug('shutting down master device')
        self.do_main_loop = False
        socket = self.context.socket(zmq.PUB)
        socket.connect(self._main_wakeup_endpoint)
        socket.send('')
            
class Worker(ZMQBase):
    def __init__(self, zmq_context, ann_endpoint, propagator):
        ZMQBase.__init__(self, zmq_context)
        
        self.propagator = propagator
        self.blocksize = 1
        
        self.ann_endpoint = ann_endpoint
        
        self.do_listen = True
        self._listen_wakeup_endpoint = 'inproc://_listen_wakeup_{:x}'.format(id(self))
        self.listen_thread = None
        
    def shutdown(self):
        log.debug('worker shutting down')
        self.do_listen = False
        socket = self.context.socket(zmq.PUB)
        socket.connect(self._listen_wakeup_endpoint)
        socket.send('')
        
    def listen(self):
        if self.listen_thread is None:
            self.listen_thread = threading.Thread(target=self.listen_loop)
            self.listen_thread.daemon = False
            self.listen_thread.start()
            
    def listen_loop(self):
        '''Receive announcements and commands from the master'''
        
        log.debug('listening')
        
        ann_socket = self.context.socket(zmq.SUB)
        ann_socket.setsockopt(zmq.SUBSCRIBE,'')
        ann_socket.connect(self.ann_endpoint)
        
        wakeup_socket = self.context.socket(zmq.SUB)
        wakeup_socket.setsockopt(zmq.SUBSCRIBE,'')
        wakeup_socket.bind(self._listen_wakeup_endpoint)

        poller = zmq.Poller()
        poller.register(ann_socket, zmq.POLLIN)
        poller.register(wakeup_socket, zmq.POLLIN)
        
        while self.do_listen:
            poll_results = set(fd for (fd,_flag) in poller.poll())
            
            if not poll_results:
                log.debug('timeout')
            
            if wakeup_socket in poll_results:
                log.debug('wake up!')
                wakeup_socket.recv()
                #self.flush_socket(wakeup_socket)
            
            if ann_socket in poll_results:
                log.debug('receiving from ann_socket')
                msg = ann_socket.recv()
                log.debug('ann_socket -> {!r}'.format(msg))
                if msg == 'shutdown':
                    log.debug('shutdown received')
                    #self.flush_socket(ann_socket)
                    return
                elif msg.startswith('tasks available '):
                    task_endpoint = msg[16:].strip()
                    tasks = True
                    
                    while tasks:
                        task_socket = self.context.socket(zmq.REQ)
                        task_socket.setsockopt(zmq.HWM,1)
                        task_socket.connect(task_endpoint)
                        task_socket.send('request {}'.format(self.blocksize))
                        results_endpoint, tasks = task_socket.recv_pyobj()
                        log.debug('received {} tasks with results destined for {}'.format(len(tasks), results_endpoint))
                        if tasks:
                            results_socket = self.context.socket(zmq.PUSH)
                            results_socket.setsockopt(zmq.HWM, 1)
                            results_socket.connect(results_endpoint)
                            for task in tasks:
                                log.debug('processing task {}'.format(task.task_id))
                                if task.task_type == 'propagate':
                                    segments = list(task.operand)
                                    self.propagator.propagate(segments)
                                    result = segments
                                log.debug('sending results for {} to {}'.format(task.task_id, results_endpoint))
                                results_socket.send_pyobj((task.task_id, result))
                                log.debug('results sent')
                        # closing and/or deleting results_socket and task_socket cause deadlocks or data loss
                else:
                    log.debug('unknown message received')

class Node(ZMQBase):
    '''Forward traffic between an upstream master (or node) and downstream workers (or nodes)'''
    def __init__(self, context,
                 upstream_task_endpoint, upstream_results_endpoint, upstream_ann_endpoint,
                 downstream_task_endpoint, downstream_results_endpoint, downstream_ann_endpoint):
        ZMQBase.__init__(self, context)
        
        # upstream (remote) endpoints
        self.upstream_task_endpoint = upstream_task_endpoint
        self.upstream_results_endpoint = upstream_results_endpoint
        self.upstream_ann_endpoint = upstream_ann_endpoint
        
        # where downstream (local) workers reach *this* object
        self.downstream_task_endpoint = downstream_task_endpoint
        self.downstream_results_endpoint = downstream_results_endpoint
        self.downstream_ann_endpoint = downstream_ann_endpoint
        
        # Start a streamer to forward results back upstream (to the master)
        self.zstreamer = zmq.devices.ThreadDevice(zmq.STREAMER, zmq.PULL, zmq.PUSH)
        self.zstreamer.bind_in(self.downstream_results_endpoint)
        self.zstreamer.connect_out(self.upstream_results_endpoint)
        log.debug('starting STREAMER')
        self.zstreamer.start()
        
        # Start a queue to forward results between upstream (the master) and downstream (clients)
        self.zqueue = zmq.devices.ThreadDevice(zmq.QUEUE, zmq.REP, zmq.REQ)
        self.zqueue.bind_in(self.downstream_task_endpoint)
        self.zqueue.connect_out(self.upstream_task_endpoint)
        log.debug('starting QUEUE')
        self.zqueue.start()
        
        # We do manual forwarding of PUB/SUB messages to catch the "shutdown" message ourselves
        
        self.forward_thread = None
        
    def forward_loop(self):
        upstream_ann = self.context.socket(zmq.SUB)
        upstream_ann.setsockopt(zmq.SUBSCRIBE, '')
        upstream_ann.connect(self.upstream_ann_endpoint)
        
        downstream_ann = self.context.socket(zmq.PUB)
        downstream_ann.bind(self.downstream_ann_endpoint)
        
        poller = zmq.Poller()
        poller.register(upstream_ann, zmq.POLLIN)
        while True:
            poll_results = set(fd for (fd, _flag) in poller.poll())
            
            if upstream_ann in poll_results:
                msg = upstream_ann.recv()
                # Immediately dispatch message
                downstream_ann.send(msg)
                
                if msg == 'shutdown':
                    return
                
    def start_forward_loop(self):
        if self.forward_thread is None:
            self.forward_thread = threading.Thread(target = self.forward_loop)
            self.forward_thread.daemon = False
            self.forward_thread.start()
    
    
                
class ZMQWorkManager(ZMQBase, WEMDWorkManager):
    def __init__(self, sim_manager):
        
        ZMQBase.__init__(self,context=None)
        WEMDWorkManager.__init__(self,sim_manager)
        
        runtime_config = sim_manager.runtime_config
        self.upstream_host  = runtime_config.get('zmq_work_manager.master', socket.gethostname())
        self.ann_port  = runtime_config.get_int('zmq_work_manager.ann_port', DEFAULT_ANN_PORT)
        self.task_port    = runtime_config.get_int('zmq_work_manager.task_port', DEFAULT_TASK_PORT)
        self.results_port = runtime_config.get_int('zmq_work_manager.results_port', DEFAULT_RESULTS_PORT)
        self.n_workers = 1
        self.blocksize = 1
        
        # Where to look for tasks, announcements, and results
        self.local_ann_endpoint = None
        self.local_task_endpoint = None
        self.local_results_endpoint = None
        
        self.remote_ann_endpoint = None
        self.remote_task_endpoint = None
        self.remote_results_endpoint = None
        
        # The ZMQ widget (Master, Worker, or Node) associated with this
        # process
        self.zdev = None

    def parse_aux_args(self, aux_args, do_help = False):
        parser = argparse.ArgumentParser(usage='%(prog)s [NON_WORK_MANAGER_OPTIONS] [OPTIONS]',
                                         add_help=False)
        
        runtime_config = self.sim_manager.runtime_config
        parser.add_argument('-n', '-np', '-nw', type=int, dest='n_workers', default=multiprocessing.cpu_count(),
                            help='Number of worker processes to run on this host. Use 0 for a dedicated server '
                                +' or forwarder process. (Default: %(default)s)')
        parser.add_argument('-H', '--host', default=socket.gethostname(),
                            help='Upstream (master) host (default: %(default)s)')
        parser.add_argument('--aport', type=int, default=DEFAULT_ANN_PORT,
                            help='Port which the master uses to send commands and announcements to '
                                +'remote workers (default: %(default)s)')
        parser.add_argument('--tport', type=int, default=DEFAULT_TASK_PORT,
                            help='Port where remote workers contact the master to receive work (default: %(default)s)')
        parser.add_argument('--rport', type=int, default=DEFAULT_RESULTS_PORT,
                            help='Port where remote workers return completed tasks to the master (default: %(default)s)')
        
        master_opts = parser.add_argument_group('options for WEMD master processes')
        master_opts.add_argument('-b', '--blocksize', type=int, dest='blocksize', default=1,
                            help='Number of segments to dispatch as a unit to each worker (default: %(default)s)')
        master_opts.add_argument('--announce-interval', type=float, dest='announce_interval', default=60,
                            help='How often (in seconds) to announce to workers that work remains to be done. '
                                +'Too low a value will increase network load. (Default: %(default)s).')
        master_opts.add_argument('--results-check-interval', type=float, dest='results_check_interval', default=0.1,
                            help='How often (in seconds) to check for results. Too low a value will increase CPU '
                                +'load on the master. (Default: %(default)s).')

        if do_help:
            parser.print_help()
            sys.exit(0)
        args, extra_args = parser.parse_known_args(aux_args)
        
        # zmq expects an IP address, not a hostname, for TCP transports
        # otherwise, the cryptic "no such device" error happens
        self.upstream_host =  runtime_config['zmq_work_manager.upstream_host'] = socket.gethostbyname(args.host)
        self.ann_port =       runtime_config['zmq_work_manager.ann_port'] = args.aport
        self.task_port =      runtime_config['zmq_work_manager.task_port'] = args.tport
        self.results_port =   runtime_config['zmq_work_manager.results_port'] = args.rport
        self.n_workers =      runtime_config['zmq_work_manager.n_workers'] = args.n_workers
        self.blocksize =      runtime_config['zmq_work_manager.blocksize'] = args.blocksize
        
        return extra_args
        
    def shutdown(self, exit_code = 0):
        log.debug('shutting down with code {}'.format(exit_code))
        if self.mode == 'master':
            self.zdev.shutdown()
        self.remove_ipc_endpoints()
                        
    def run_worker(self):
        assert self.mode in ('worker', 'nodeworker')
        assert self.zdev.listen_thread is not None
        return self.zdev.listen_thread.join()
    
    def run_node(self):
        assert self.mode == 'node'
        assert self.zdev.forward_thread is not None
        return self.zdev.forward_thread.join()
    
    def spawn_worker_processes(self):
        if self.mode == 'master':
            n_to_spawn = self.n_workers
        elif self.mode == 'worker':
            # We started as a dedicated worker, so this process *is* a worker process
            # and we only need to spawn enough processes to fill up to self.n_workers
            n_to_spawn = self.n_workers - 1
        elif self.mode == 'node':
            n_to_spawn = self.n_workers
        else:
            raise ValueError('invalid operational mode {!r}'.format(self.mode))
        
        if n_to_spawn <= 0:
            return
        
        for iworker in xrange(0, n_to_spawn):
            log.debug('spawning worker {} of {}'.format(iworker+1, n_to_spawn))
            pid = os.fork()
            if pid == 0:
                # in child
                if self.mode == 'master':
                    self.mode = 'worker'
                elif self.mode == 'node':
                    self.mode = 'nodeworker'
                return
            
    def prepare(self):
        
        # Assign paths for local IPC endpoints; these will be carried across fork()s
        self.local_ann_endpoint = self.make_ipc_endpoint()
        self.local_task_endpoint = self.make_ipc_endpoint()
        self.local_results_endpoint = self.make_ipc_endpoint()
        
        # Connect to a specific remote host for now; maybe add wildcards for the master later
        self.remote_ann_endpoint = 'tcp://{}:{}'.format(self.upstream_host, self.ann_port)
        self.remote_task_endpoint = 'tcp://{}:{}'.format(self.upstream_host, self.task_port)
        self.remote_results_endpoint = 'tcp://{}:{}'.format(self.upstream_host, self.results_port)
                
        # Spawn workers; we return here after fork()s, with self.mode reset appropriately
        # (i.e. worker processes have self.mode set to 'worker')
        self.spawn_worker_processes()
        
        assert self.context is None
        assert self.zdev is None
        
        # Create ZMQ context for this process
        self.context = zmq.Context()
        
        # And configure and run the ZMQ widget for this process
        if self.mode == 'master':
            log.debug('PID {} is master'.format(os.getpid()))            
            self.zdev = Master(self.context, self.remote_task_endpoint, self.remote_results_endpoint, self.remote_ann_endpoint)
            #self.zdev = Master(self.context, self.local_task_endpoint, self.local_results_endpoint, self.local_ann_endpoint)
            self.zdev.start_main_thread()
        elif self.mode == 'node':
            log.debug('PID {} is node coordinator'.format(os.getpid()))
            self.zdev = Node(self.context, self.remote_task_endpoint, self.remote_results_endpoint, self.remote_ann_endpoint,
                                           self.local_task_endpoint, self.local_results_endpoint, self.local_ann_endpoint)
            self.zdev.start_forward_loop()
        elif self.mode == 'worker':
            log.debug('PID {} is worker'.format(os.getpid()))
            self.zdev = Worker(self.context, self.remote_ann_endpoint, self.sim_manager.propagator)
            #self.zdev = Worker(self.context, self.local_ann_endpoint, self.sim_manager.propagator)
            # Listen for work immediately; run_worker() then joins the worker thread to await termination
            self.zdev.listen()
        elif self.mode == 'nodeworker':
            log.debug('PID {} is node worker'.format(os.getpid()))
            self.zdev = Worker(self.context, self.local_ann_endpoint, self.sim_manager.propagator)
            self.zdev.listen()
        else:
            raise NotImplementedError
            
    
    def propagate(self, segments):
        assert self.mode == 'master'
        
        blocks = [segments[i:i+self.blocksize] for i in xrange(0,len(segments),self.blocksize)]
        log.debug('dispatching {} segment(s) in {} block(s)'.format(len(segments), len(blocks)))
        
        # label deterministically to save a little system entropy
        tasks = [Task('propagate', block, task_id=('propagate', block[0].n_iter, block[0].seg_id)) for block in blocks]
        self.zdev.dispatch_all(tasks)
        results = self.zdev.wait_all(task.task_id for task in tasks)
        
        # No unprocessed results?
        assert len(self.zdev.results) == 0
        
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
