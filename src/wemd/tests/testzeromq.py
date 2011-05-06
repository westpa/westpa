from __future__ import division, print_function; __metaclass__ = type
import threading, time
import cPickle as pickle
import wemd.work_managers.zeromq
import nose.tools
import zmq

from wemd.work_managers.zeromq import ZWorkProvider


TIME_FAST = 0.5
TIME_MEDIUM = 4*TIME_FAST
TIME_SLOW = 5*TIME_MEDIUM

DELAY_FAST = TIME_FAST/10.0
DELAY_MEDIUM = TIME_MEDIUM/10.0
DELAY_SLOW = TIME_SLOW/10.0


class TestZWorkProvider:
    def __init__(self):
        self.context = zmq.Context()
        
        self.task_endpoint = 'inproc://tasks'
        self.results_endpoint = 'inproc://results'
        self.ann_endpoint = 'inproc://announce'
        self.ann_socket = self.context.socket(zmq.SUB)
        self.ann_socket.setsockopt(zmq.SUBSCRIBE,'')
        self.ann_socket.bind(self.ann_endpoint)
        
    def make_provider(self):
        zwp = ZWorkProvider(self.context, self.task_endpoint, self.results_endpoint, self.ann_endpoint)
        return zwp
            
    @nose.tools.timed(TIME_FAST)
    def test_announce(self):        
        zwp = self.make_provider()
        zwp.announce_work()
        rsl = self.ann_socket.recv()
        assert rsl.startswith('tasks available')

        
    @nose.tools.timed(TIME_FAST)
    def test_dispatch(self):
        zwp = self.make_provider()
        
        task_socket = self.context.socket(zmq.REQ)
        task_socket.connect(self.task_endpoint)
        
        #zwp.outgoing_tasks.extend([1,'2', 3.0])

        #t = threading.Thread(target=zwp.dispatch_loop)
        #t.start()
        zwp.dispatch([1,'2', 3.0])
        rsl = self.ann_socket.recv()
        
        task_socket.send('request 2')
        rsl = task_socket.recv_multipart()
        assert rsl[0] == 'tasks 2 {}'.format(self.results_endpoint)
        payload = pickle.loads(rsl[1]) 
        assert payload == [1,'2']
        
        task_socket.send('request 2')
        rsl = task_socket.recv_multipart()
        assert rsl[0] == 'tasks 1 {}'.format(self.results_endpoint)
        payload = pickle.loads(rsl[1])
        assert payload == [3.0]
        
        zwp.dispatch([])
        task_socket.send('request 1')
        rsl = task_socket.recv_multipart()
        assert rsl[0] == 'tasks 0'
        assert rsl[1] == ''
        zwp.shutdown()
        
    @nose.tools.timed(TIME_FAST)
    def test_collect(self):
        zwp = self.make_provider()
        zwp.collect()
        
        results_socket = self.context.socket(zmq.PUSH)
        results_socket.connect(self.results_endpoint)
        
        results_socket.send('propagate', zmq.SNDMORE)
        results_socket.send_pyobj([1, '2', 3.0])
        time.sleep(DELAY_FAST)
        assert list(zwp.incoming_results) == [1, '2', 3.0]
        zwp.shutdown()
        