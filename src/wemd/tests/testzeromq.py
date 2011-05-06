from __future__ import division, print_function; __metaclass__ = type
import threading, time
import cPickle as pickle
import wemd.work_managers.zeromq
import nose.tools
import zmq

from wemd.work_managers.zeromq import ZWorkProvider, ZWorkConsumer


TIME_FAST = 0.5
TIME_MEDIUM = 2
TIME_SLOW = 10

DELAY_FAST = 0.05
DELAY_MEDIUM = 0.2
DELAY_SLOW = 1


class TestZWorkProvider:
    def __init__(self):
        self.context = zmq.Context()
        
        #self.task_endpoint = 'inproc://tasks'
        #self.results_endpoint = 'inproc://results'
        #self.ann_endpoint = 'inproc://announce'
        
        #self.task_endpoint = 'ipc:///tmp/tasks'
        #self.results_endpoint = 'ipc:///tmp/results'
        #self.ann_endpoint = 'ipc:///tmp/announce'
        
        self.task_endpoint = 'tcp://127.0.0.1:23811'
        self.results_endpoint = 'tcp://127.0.0.1:23812'
        self.ann_endpoint = 'tcp://127.0.0.1:23813'
                
    def make_provider(self):
        zwp = ZWorkProvider(self.context, self.task_endpoint, self.results_endpoint, self.ann_endpoint)
        self.ann_socket = self.context.socket(zmq.SUB)
        self.ann_socket.setsockopt(zmq.SUBSCRIBE,'')
        self.ann_socket.setsockopt(zmq.LINGER,0)
        self.ann_socket.connect(self.ann_endpoint)
        time.sleep(DELAY_FAST)
        return zwp
            
    @nose.tools.timed(TIME_FAST)
    def test_announce(self):    
        zwp = self.make_provider()
        zwp.announce_work()
        rsl = self.ann_socket.recv()
        assert rsl.startswith('tasks available')
        
        # A delay after shutdown appears necessary before running
        # other tests
        zwp.shutdown()
        time.sleep(DELAY_FAST)

        
    @nose.tools.timed(TIME_MEDIUM)
    def test_dispatch(self):
        zwp = self.make_provider()
        
        task_socket = self.context.socket(zmq.REQ)
        task_socket.connect(self.task_endpoint)

        zwp.dispatch([1,'2', 3.0])
        
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
        time.sleep(DELAY_FAST)
                
    @nose.tools.timed(TIME_FAST)
    def test_collect(self):
        zwp = self.make_provider()
        zwp.collect()
        
        results_socket = self.context.socket(zmq.PUSH)
        results_socket.connect(self.results_endpoint)
        results_socket.send_pyobj([1, '2', 3.0])
        time.sleep(DELAY_FAST)
        assert list(zwp.incoming_results) == [1, '2', 3.0]

        zwp.shutdown()
        time.sleep(DELAY_FAST)

class TestZWorkConsumer:
    def __init__(self):
        self.context = zmq.Context()
        
        self.task_endpoint = 'inproc://tasks'
        self.results_endpoint = 'inproc://results'
        self.ann_endpoint = 'inproc://announce'
        
    def make_provider(self):
        zwp = ZWorkProvider(self.context, self.task_endpoint, self.results_endpoint, self.ann_endpoint)
        time.sleep(DELAY_FAST)
        return zwp
    
    def make_consumer(self, blocksize, propagator):
        zwc = ZWorkConsumer(self.context, self.ann_endpoint, blocksize, propagator)
        time.sleep(DELAY_FAST)
        return zwc
    
    def test_work(self):
        class Propagator:
            def propagate(self, segments):
                for i in xrange(0, len(segments)):
                    segment = segments[i]
                    segment = [pcoord + 1.0 for pcoord in segment]
                    segments[i] = segment
        
        zwp = self.make_provider()
        zwp.dispatch([('propagate', [1.0,2.0,3.0-1j]), ('propagate', [5])])
        zwp.collect()
        
        zwc = self.make_consumer(2, Propagator())
        zwc.work()
        nloops = 0
        while nloops < 10 and len(zwp.incoming_results) == 0:
            time.sleep(DELAY_FAST)
        assert nloops != 10    
        assert list(zwp.incoming_results) == [('propagate', [2.0, 3.0, 4.0-1j]),
                                              ('propagate', [6])]
        
        zwc.shutdown()
        zwp.shutdown()            