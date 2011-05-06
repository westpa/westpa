from __future__ import division, print_function; __metaclass__ = type

import sys, os, tempfile, logging, socket, multiprocessing, cPickle, threading, time, uuid, signal
import collections
import zmq

log = logging.getLogger(__name__)

class ZMQBase:
    def __init__(self, zmq_context):
        self.context = zmq_context
        self.hostname = socket.gethostname()
        self.nodespec = '{}#{:d}:0x{:x}'.format(self.hostname, os.getpid(), threading.current_thread().ident)
        self.node_id = str(uuid.uuid4())        
        
        # An inproc PUB socket for asynchronous inter-thread signaling
        self._intsig_endpoint = None
        self._intsig = None
        
    def bind_intsig(self):
        self._intsig_endpoint = 'inproc://_intsig'
        self._intsig = self.context.socket(zmq.PUB)
        self._intsig.bind(self._intsig_endpoint)
        
    def listen_intsig(self, prefix=''):
        if self._intsig is None:
            self.bind_intsig()
        socket = self.context.socket(zmq.SUB)
        socket.setsockopt(zmq.SUBSCRIBE, prefix)
        socket.connect(self._intsig_endpoint)
        return socket
    
    def intsig(self, message):
        '''Send an in-process signal to all registered listeners'''
        if self._intsig is None:
            # No one listening, but bind anyway
            self.bind_intsig()
        self._intsig.send(message)
    
        
class ZWorkProvider(ZMQBase):
    def __init__(self, zmq_context, task_endpoint, results_endpoint, ann_out_endpoint):
        ZMQBase.__init__(self, zmq_context)
        
        self.outgoing_tasks = collections.deque()
        self.incoming_results = collections.deque()
        
        self.task_endpoint = task_endpoint
        self.results_endpoint = results_endpoint
        self.ann_out_endpoint = ann_out_endpoint
        
        self.task_socket = self.context.socket(zmq.REP)
        self.results_socket = self.context.socket(zmq.PULL)
        self.ann_out_socket = None
        
        self.task_socket.bind(task_endpoint)
        self.results_socket.bind(results_endpoint)
        
        self.do_dispatch = True
        self.dispatch_thread = None
        
        self.do_collect = True
        self.collect_thread = None
        
        # Timeouts are fatal
        self.test_mode = False
                
    def announce_work(self):
        if self.ann_out_socket is None:
            self.ann_out_socket = self.context.socket(zmq.PUB)
            self.ann_out_socket.connect(self.ann_out_endpoint)
        self.ann_out_socket.send('tasks available {}'.format(self.task_endpoint))
        
    def handle_work_request(self):
        req = self.task_socket.recv()
        if not req.startswith('request'):
            return
        
        _reqtext, ntasktext = req.split()
        n_tasks = int(ntasktext)
        to_send = []
        while len(to_send) < n_tasks:
            try:
                to_send.append(self.outgoing_tasks.popleft())
            except IndexError:
                break        
        if to_send:
            self.task_socket.send('tasks {} {}'.format(len(to_send), self.results_endpoint), zmq.SNDMORE)
            self.task_socket.send_pyobj(to_send)
        else:
            self.task_socket.send('tasks 0', zmq.SNDMORE)
            self.task_socket.send('')
            
    def handle_results(self):
        rsl = self.results_socket.recv_multipart()
        if rsl[0] not in ('propagate',):
            return
        payload = cPickle.loads(rsl[1])
        for result in payload:
            self.incoming_results.append(result)
        
    def dispatch(self, tasks):
        for task in tasks:
            self.outgoing_tasks.append(task)
        if self.outgoing_tasks:
            self.announce_work()
            
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
        poller = zmq.Poller()
        poller.register(self.task_socket, zmq.POLLIN)
        poller.register(intsig, zmq.POLLIN)
        
        if self.outgoing_tasks:
            self.announce_work()
            
        while self.do_dispatch:
            poll_results = set(poller.poll())
            if (self.task_socket, zmq.POLLIN) in poll_results:
                self.handle_work_request()
            if (intsig, zmq.POLLIN) in poll_results and intsig.recv() == 'exit':
                return
            
    def collect_loop(self):
        intsig = self.listen_intsig()
        poller = zmq.Poller()
        poller.register(self.results_socket, zmq.POLLIN)
        poller.register(intsig, zmq.POLLIN)
        
        while self.do_collect:
            poll_results = set(poller.poll())
            if (self.results_socket, zmq.POLLIN) in poll_results:
                self.handle_results()
            if (intsig, zmq.POLLIN) in poll_results and intsig.recv() == 'exit':
                return
            
class ZWorkConsumer(ZMQBase):
    def __init__(self, zmq_context, ann_in_endpoint,
                 blocksize, propagator):
        ZMQBase.__init__(self, zmq_context)
        
        self.blocksize = blocksize
        self.propagator = propagator
        
        self.ann_in_endpoint = ann_in_endpoint
        
        self.ann_in_socket = self.context.bind(zmq.SUB)
        self.ann_in_socket.setsockopt(zmq.SUBSCRIBE,'')
        self.ann_in_socket.bind(self.ann_in_endpoint)
                
        self.do_work = True
    
    def work_loop(self):
        intsig = self.listen_intsig()
        poller = zmq.Poller()
        poller.register(self.ann_in_socket, zmq.POLLIN)
        poller.register(intsig, zmq.POLLIN)
        
        while self.do_work:
            poll_results = set(poller.poll())
            
            if (intsig, zmq.POLLIN) in poll_results and intsig.recv() == 'exit':
                return
            
            if (self.ann_in_socket, zmq.POLLIN) in poll_results:
                msg = self.ann_in_socket.recv()
                if not msg.startswith('tasks available'):
                    log.error('invalid reply to task available query; ignoring') 
                    continue
                _prefix, task_endpoint = msg.rsplit(None, 1)
                
                task_socket = self.context.connect(task_endpoint)
                task_socket.send('request {}'.format(self.blocksize))
                rpl = task_socket.recv_multipart()
                if not rpl[0].startswith('tasks'):
                    log.error('invalid reply to task request; ignoring')
                    continue
                
                _prefix, ntaskstr, results_endpoint = rpl[0].split()
                n_tasks = int(ntaskstr)
                if not n_tasks: 
                    continue
                elif not (n_tasks <= self.blocksize):
                    log.error('invalid number of tasks received; ignoring')
                    continue
                
                payload = cPickle.loads(rpl[1])
                if len(payload) != n_tasks:
                    log.error('payload contained wrong number of tasks; ignoring')
                    continue
                
                
                
                
                
                
        