# Copyright (C) 2013 Matthew C. Zwier, Nick Rego, and Lillian T. Chong
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
import os, signal, tempfile, time, sys, multiprocessing, uuid, socket, json, re
import cPickle as pickle

from work_managers import WMFuture
from work_managers import zeromq as zwm
from work_managers.zeromq import ZMQWorkManager, ZMQClient, recvall, ZMQServer, ZMQRouter
from work_managers.zeromq.client import ZMQWMProcess
from work_managers.zeromq.router import ZMQDevice
from tsupport import *

import zmq

import nose.tools
from nose.tools import raises, nottest, timed, assert_raises #@UnresolvedImport
from nose.plugins.skip import SkipTest

def sockdelay():
    '''Delay for slightly longer than the default auto-reconnect time for ZeroMQ (100 ms)'''
    time.sleep(0.11)
    
def randport():
    s = socket.socket()
    s.bind(('127.0.0.1',0))
    port = s.getsockname()[1]
    s.close()
    return port

def randipc():
    (fd, socket_path) = tempfile.mkstemp()
    os.close(fd)
    endpoint = 'ipc://{}'.format(socket_path)
    return endpoint

class BaseTestZMQEnvironment:
    sanitize_vars = ('WM_WORK_MANAGER', 'WM_N_WORKERS',
                     'WM_ZMQ_INFO','WM_ZMQ_WRITE_INFO', 'WM_ZMQ_READ_INFO', 'WM_ZMQ_CLIENT_COMM_MODE',
                     'WM_ZMQ_TASK_ENDPOINT', 'WM_ZMQ_RESULT_ENDPOINT', 'WM_ZMQ_ANNOUNCE_ENDPOINT', 'WM_ZMQ_LISTEN_ENDPOINT',
                     'WM_ZMQ_DOWNSTREAM_TASK_ENDPOINT', 'WM_ZMQ_DOWNSTREAM_RESULT_ENDPOINT', 'WM_ZMQ_DOWNSTREAM_ANNOUNCE_ENDPOINT',
                     'WM_ZMQ_UPSTREAM_TASK_ENDPOINT', 'WM_ZMQ_UPSTREAM_RESULT_ENDPOINT', 'WM_ZMQ_UPSTREAM_ANNOUNCE_ENDPOINT',
                     'WM_ZMQ_UPSTREAM_LISTEN_ENDPOINT', 'WM_ZMQ_DOWNSTREAM_LISTEN_ENDPOINT')
    
    def setUp(self):
        for varname in self.sanitize_vars:
            assert varname not in os.environ
 
    def tearDown(self):
        for varname in self.sanitize_vars:
            os.environ.pop(varname, None)

class TestZMQRouterEnvironment(BaseTestZMQEnvironment):

    def setUp(self):
        super(TestZMQRouterEnvironment, self).setUp()
        #assert 'WM_ZMQ_WRITE_INFO' not in os.environ
        #assert 'WM_ZMQ_DOWNSTREAM_TASK_ENDPOINT' not in os.environ
        self.filename = 'testfile'
        self.router = None

        self.up_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.up_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.up_announce_endpoint = 'tcp://127.0.0.1:{}'.format(randport())  
        self.up_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        
        self.hostname = socket.gethostname()
        with open(self.filename, 'wt') as infofile:
            json.dump({'task_endpoint': re.sub(r'\*', self.hostname, self.up_task_endpoint),
                       'result_endpoint': re.sub(r'\*', self.hostname, self.up_result_endpoint),
                       'announce_endpoint': re.sub(r'\*', self.hostname, self.up_announce_endpoint),
                       'listen_endpoint': re.sub(r'\*', self.hostname, self.up_listen_endpoint)},
                      infofile)
        os.chmod(self.filename, 0600)

    def tearDown(self):
        super(TestZMQRouterEnvironment, self).tearDown()
        try:
            os.unlink(self.filename)
        finally:
            if self.router is not None:
                self.router.remove_router_info_file()
                del self.router

    def test_read_upstream_info(self):
        '''Router environment: Reads from info file correctly when no upstream endpoints specified'''
        os.environ['WM_ZMQ_READ_INFO'] = self.filename
        self.router = ZMQRouter.from_environ()

        assert self.router.upstream_task_endpoint == self.up_task_endpoint, 'router upstream task endpoint does not match server file'
        assert self.router.upstream_result_endpoint == self.up_result_endpoint, 'router upstream result endpoint does not match server file'
        assert self.router.upstream_announce_endpoint == self.up_announce_endpoint, 'router upstream announce endpoint does not match file'

    def test_write_router_info(self):
        '''Router environment: Writes router info file successfully'''

        os.environ['WM_ZMQ_READ_INFO'] = self.filename

        self.router = ZMQRouter.from_environ()

        assert self.router.router_info_filename is not None

        router_info = json.load(open(self.router.router_info_filename,'rt'))

        assert re.sub(r'\*', self.hostname, self.router.downstream_task_endpoint) == router_info['task_endpoint']
        assert re.sub(r'\*', self.hostname, self.router.downstream_result_endpoint) == router_info['result_endpoint']
        assert re.sub(r'\*', self.hostname, self.router.downstream_announce_endpoint) == router_info['announce_endpoint']

    @raises(EnvironmentError)
    def test_router_bad_upstream_environment(self):
        '''Router environment: Raises EnvironmentError if upstream task endpoints and server info file not specified'''

        self.router = ZMQRouter.from_environ()

class TestZMQWorkManager(BaseTestZMQEnvironment):

    def test_auto_local(self):
        with ZMQWorkManager() as work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()

    def test_environ_empty(self):
        with ZMQWorkManager.from_environ() as work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()

    def test_server_info_ipc(self):
        with ZMQWorkManager() as work_manager:
            server_info_filename = work_manager.server_info_filename
            assert os.path.exists(server_info_filename)
            assert os.stat(server_info_filename).st_mode & 0777 == 0600
            server_info = json.load(open(server_info_filename, 'rt'))
            assert re.sub(r'\*', socket.gethostname(), work_manager.master_task_endpoint) == server_info['task_endpoint']
            assert re.sub(r'\*', socket.gethostname(), work_manager.master_result_endpoint) == server_info['result_endpoint']
            assert re.sub(r'\*', socket.gethostname(), work_manager.master_announce_endpoint) == server_info['announce_endpoint']
            assert re.sub(r'\*', socket.gethostname(), work_manager.master_listen_endpoint) == server_info['listen_endpoint']
        assert not os.path.exists(server_info_filename)

    def test_environ_nworkers(self):
        os.environ['WM_N_WORKERS'] = str(2)
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.internal_client.n_workers == 2
            future = work_manager.submit(will_succeed)
            future.get_result()

    def test_environ_noworkers(self):
        os.environ['WM_N_WORKERS'] = str(0)
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.internal_client is None
                    
    def test_worker_ids(self):
        work_manager = ZMQWorkManager()
        with work_manager:
            futures = work_manager.submit_many([(get_process_index, (), {})] * work_manager.n_local_workers)
            work_manager.wait_all(futures)
            results = set(future.get_result() for future in futures)
            assert results == set(str(n) for n in xrange(work_manager.n_local_workers))
            
    @raises(ValueError)
    def test_client_from_bad_environ(self):
        os.environ['WM_N_WORKERS'] = str(0)
        task_endpoint = 'tcp://*:{}'.format(randport())
        result_endpoint = 'tcp://*:{}'.format(randport())
        announce_endpoint = 'tcp://*:{}'.format(randport())
        listen_endpoint = 'tcp://*:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        os.environ['WM_ZMQ_LISTEN_ENDPOINT'] = listen_endpoint
        
        with ZMQWorkManager() as work_manager:
            os.environ['WM_N_WORKERS'] = str(2)
            test_client = ZMQClient.from_environ()
            test_client.startup()            
            try:
                future = work_manager.submit(will_succeed)
                future.get_result()                
            finally:
                test_client.shutdown()

    def test_client_from_environ(self):
        os.environ['WM_N_WORKERS'] = str(0)
        task_endpoint = 'tcp://localhost:{}'.format(randport())
        result_endpoint = 'tcp://localhost:{}'.format(randport())
        announce_endpoint = 'tcp://localhost:{}'.format(randport())
        listen_endpoint = 'tcp://localhost:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        os.environ['WM_ZMQ_LISTEN_ENDPOINT'] = listen_endpoint
        
        with ZMQWorkManager() as work_manager:
            os.environ['WM_N_WORKERS'] = str(2)
            test_client = ZMQClient.from_environ()
            test_client.startup()            
            try:
                future = work_manager.submit(will_succeed)
                future.get_result()                
            finally:
                test_client.shutdown()
                
    def test_client_from_environ_tcp(self):
        os.environ['WM_N_WORKERS'] = str(0)
        task_endpoint = 'tcp://localhost:{}'.format(randport())
        result_endpoint = 'tcp://localhost:{}'.format(randport())
        announce_endpoint = 'tcp://localhost:{}'.format(randport())
        listen_endpoint = 'tcp://localhost:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        os.environ['WM_ZMQ_LISTEN_ENDPOINT'] = listen_endpoint
        os.environ['WM_ZMQ_CLIENT_COMM_MODE'] = 'tcp'
        
        with ZMQWorkManager() as work_manager:
            os.environ['WM_N_WORKERS'] = str(2)
            test_client = ZMQClient.from_environ()
            test_client.startup()
            assert test_client.worker_task_endpoint.startswith('tcp://')
            assert test_client.worker_result_endpoint.startswith('tcp://')
            try:
                future = work_manager.submit(will_succeed)
                future.get_result()                
            finally:
                test_client.shutdown()
                
    def test_client_from_environ_ipc(self):
        os.environ['WM_N_WORKERS'] = str(0)
        task_endpoint = 'tcp://localhost:{}'.format(randport())
        result_endpoint = 'tcp://localhost:{}'.format(randport())
        announce_endpoint = 'tcp://localhost:{}'.format(randport())
        listen_endpoint = 'tcp://localhost:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        os.environ['WM_ZMQ_LISTEN_ENDPOINT'] = listen_endpoint
        os.environ['WM_ZMQ_CLIENT_COMM_MODE'] = 'ipc'
        
        with ZMQWorkManager() as work_manager:
            os.environ['WM_N_WORKERS'] = str(2)
            test_client = ZMQClient.from_environ()
            test_client.startup()
            assert test_client.worker_task_endpoint.startswith('ipc://')
            assert test_client.worker_result_endpoint.startswith('ipc://')
            try:
                future = work_manager.submit(will_succeed)
                future.get_result()                
            finally:
                test_client.shutdown()

            
    def test_environ_tcp_endpoints(self):
    
        # note that this tests not only that the work manager honor our environment settings, but that
        # the hostname-to-ip mapping succeeded
        task_endpoint = 'tcp://localhost:{}'.format(randport())
        result_endpoint = 'tcp://localhost:{}'.format(randport())
        announce_endpoint = 'tcp://localhost:{}'.format(randport())
        listen_endpoint = 'tcp://localhost:{}'.format(randport())

        os.environ['WM_ZMQ_TASK_ENDPOINT'] = task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = announce_endpoint
        os.environ['WM_ZMQ_LISTEN_ENDPOINT'] = listen_endpoint
        
        with ZMQWorkManager.from_environ() as work_manager:
            assert work_manager.master_task_endpoint == re.sub('localhost','127.0.0.1', task_endpoint)
            assert work_manager.master_result_endpoint == re.sub('localhost','127.0.0.1', result_endpoint)
            assert work_manager.master_announce_endpoint == re.sub('localhost','127.0.0.1', announce_endpoint)
            assert work_manager.master_listen_endpoint == re.sub('localhost','127.0.0.1', listen_endpoint)
            future = work_manager.submit(will_succeed)
            future.get_result()

class TestZMQCoordinatedEnvironNoRouter(BaseTestZMQEnvironment):

    def setUp(self):
        super(TestZMQCoordinatedEnvironNoRouter, self).setUp()
        self.master_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.master_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.master_announce_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.master_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
    
    def test_explicit_endpoint(self):
        '''Coordinated environment (no router): server and clients start up when endpoints specified in environment'''
        os.environ['WM_ZMQ_TASK_ENDPOINT'] = self.master_task_endpoint
        os.environ['WM_ZMQ_RESULT_ENDPOINT'] = self.master_result_endpoint
        os.environ['WM_ZMQ_ANNOUNCE_ENDPOINT'] = self.master_announce_endpoint
        os.environ['WM_ZMQ_LISTEN_ENDPOINT'] = self.master_listen_endpoint

        client = ZMQClient.from_environ()
        with ZMQWorkManager.from_environ() as server:

            assert client.upstream_task_endpoint == self.master_task_endpoint
            assert client.upstream_result_endpoint == self.master_result_endpoint
            assert client.upstream_announce_endpoint == self.master_announce_endpoint
            assert client.upstream_listen_endpoint == self.master_listen_endpoint

            assert server.master_task_endpoint == self.master_task_endpoint
            assert server.master_result_endpoint == self.master_result_endpoint
            assert server.master_announce_endpoint == self.master_announce_endpoint
            assert server.master_listen_endpoint == self.master_listen_endpoint
    
    def test_noexplicit_endpoint(self):
        '''Coordinated environment (no router): server and clients start up when no explicit endpoints specified (i.e. via info file)'''
        
        filename = 'server_info_test'
        hostname = socket.gethostname()
        os.environ['WM_ZMQ_INFO'] = filename

        with ZMQWorkManager.from_environ() as server:

            #Wait for file to be written
            while not os.access(filename, os.F_OK):
                time.sleep(.11)

            client = ZMQClient.from_environ()

            assert re.sub(r'\*', hostname, server.master_task_endpoint) == client.upstream_task_endpoint
            assert re.sub(r'\*', hostname, server.master_result_endpoint) == client.upstream_result_endpoint
            assert re.sub(r'\*', hostname, server.master_announce_endpoint) == client.upstream_announce_endpoint
            assert re.sub(r'\*', hostname, server.master_listen_endpoint) == client.upstream_listen_endpoint

class TestZMQCoordinatedEnvironWithRouter(BaseTestZMQEnvironment):

    def setUp(self):
        super(TestZMQCoordinatedEnvironWithRouter, self).setUp()

        self.server_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.server_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.server_announce_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.server_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())

        self.client_task_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.client_result_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.client_announce_endpoint = 'tcp://127.0.0.1:{}'.format(randport())
        self.client_listen_endpoint = 'tcp://127.0.0.1:{}'.format(randport())

    def test_all_explicit_endpoint(self):
        '''Coordinated enviro (w/ router):  successful initialization when all endpoints explicitly specified in environment'''

        os.environ['WM_ZMQ_UPSTREAM_TASK_ENDPOINT'] = self.server_task_endpoint
        os.environ['WM_ZMQ_UPSTREAM_RESULT_ENDPOINT'] = self.server_result_endpoint
        os.environ['WM_ZMQ_UPSTREAM_ANNOUNCE_ENDPOINT'] = self.server_announce_endpoint
        os.environ['WM_ZMQ_UPSTREAM_LISTEN_ENDPOINT'] = self.server_listen_endpoint

        os.environ['WM_ZMQ_DOWNSTREAM_TASK_ENDPOINT'] = self.client_task_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_RESULT_ENDPOINT'] = self.client_result_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_ANNOUNCE_ENDPOINT'] = self.client_announce_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_LISTEN_ENDPOINT'] = self.client_listen_endpoint

        ##Initialize the router from environmental variables
        router = ZMQRouter.from_environ()

        ##Note that server and client have different definitions of their environmental variables
        os.environ['WM_ZMQ_DOWNSTREAM_TASK_ENDPOINT'] = self.server_task_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_RESULT_ENDPOINT'] = self.server_result_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_ANNOUNCE_ENDPOINT'] = self.server_announce_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_LISTEN_ENDPOINT'] = self.server_listen_endpoint

        os.environ['WM_ZMQ_UPSTREAM_TASK_ENDPOINT'] = self.client_task_endpoint
        os.environ['WM_ZMQ_UPSTREAM_RESULT_ENDPOINT'] = self.client_result_endpoint
        os.environ['WM_ZMQ_UPSTREAM_ANNOUNCE_ENDPOINT'] = self.client_announce_endpoint
        os.environ['WM_ZMQ_UPSTREAM_LISTEN_ENDPOINT'] = self.client_listen_endpoint

        with ZMQWorkManager.from_environ() as server:

            client = ZMQClient.from_environ()

            #Test that upstream/downstream endpoints agree
            assert router.upstream_task_endpoint != router.downstream_task_endpoint
            assert router.upstream_result_endpoint != router.downstream_result_endpoint
            assert router.upstream_announce_endpoint != router.downstream_announce_endpoint
            assert router.upstream_listen_endpoint != router.downstream_listen_endpoint

            assert router.upstream_task_endpoint == server.master_task_endpoint == self.server_task_endpoint
            assert router.upstream_result_endpoint == server.master_result_endpoint == self.server_result_endpoint
            assert router.upstream_announce_endpoint == server.master_announce_endpoint == self.server_announce_endpoint
            assert router.upstream_listen_endpoint == server.master_listen_endpoint == self.server_listen_endpoint

            assert router.downstream_task_endpoint == client.upstream_task_endpoint == self.client_task_endpoint
            assert router.downstream_result_endpoint == client.upstream_result_endpoint == self.client_result_endpoint
            assert router.downstream_announce_endpoint == client.upstream_announce_endpoint == self.client_announce_endpoint
            assert router.downstream_listen_endpoint == client.upstream_listen_endpoint == self.client_listen_endpoint

            router.remove_router_info_file()

    def test_server_explicit_endpoint(self):
        '''Coordinated enviro (w/ router): Successful initialization when only server endpoints explicitly specified'''


        os.environ['WM_ZMQ_DOWNSTREAM_TASK_ENDPOINT'] = self.server_task_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_RESULT_ENDPOINT'] = self.server_result_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_ANNOUNCE_ENDPOINT'] = self.server_announce_endpoint
        os.environ['WM_ZMQ_DOWNSTREAM_LISTEN_ENDPOINT'] = self.server_listen_endpoint

        filename = 'router_info_test'
        hostname = socket.gethostname()

        with ZMQWorkManager.from_environ() as server:

            os.environ.pop('WM_ZMQ_DOWNSTREAM_TASK_ENDPOINT', None)
            os.environ.pop('WM_ZMQ_DOWNSTREAM_RESULT_ENDPOINT', None)
            os.environ.pop('WM_ZMQ_DOWNSTREAM_ANNOUNCE_ENDPOINT', None)
            os.environ.pop('WM_ZMQ_DOWNSTREAM_LISTEN_ENDPOINT', None)

            #client initializes from router info file
            os.environ['WM_ZMQ_WRITE_INFO'] = filename

            os.environ['WM_ZMQ_UPSTREAM_TASK_ENDPOINT'] = self.server_task_endpoint
            os.environ['WM_ZMQ_UPSTREAM_RESULT_ENDPOINT'] = self.server_result_endpoint
            os.environ['WM_ZMQ_UPSTREAM_ANNOUNCE_ENDPOINT'] = self.server_announce_endpoint
            os.environ['WM_ZMQ_UPSTREAM_LISTEN_ENDPOINT'] = self.server_listen_endpoint

            router = ZMQRouter.from_environ()

            os.environ['WM_ZMQ_READ_INFO'] = filename

            os.environ.pop('WM_ZMQ_UPSTREAM_TASK_ENDPOINT', None)
            os.environ.pop('WM_ZMQ_UPSTREAM_RESULT_ENDPOINT', None)
            os.environ.pop('WM_ZMQ_UPSTREAM_ANNOUNCE_ENDPOINT', None)
            os.environ.pop('WM_ZMQ_UPSTREAM_LISTEN_ENDPOINT', None)

            client = ZMQClient.from_environ()

            assert router.upstream_task_endpoint != router.downstream_task_endpoint
            assert router.upstream_result_endpoint != router.downstream_result_endpoint
            assert router.upstream_announce_endpoint != router.downstream_announce_endpoint
            assert router.upstream_listen_endpoint != router.downstream_listen_endpoint

            assert router.upstream_task_endpoint == server.master_task_endpoint == self.server_task_endpoint
            assert router.upstream_result_endpoint == server.master_result_endpoint == self.server_result_endpoint
            assert router.upstream_announce_endpoint == server.master_announce_endpoint == self.server_announce_endpoint
            assert router.upstream_listen_endpoint == server.master_listen_endpoint == self.server_listen_endpoint

            assert re.sub(r'\*', hostname, router.downstream_task_endpoint) == client.upstream_task_endpoint != self.client_task_endpoint
            assert re.sub(r'\*', hostname, router.downstream_result_endpoint) == client.upstream_result_endpoint != self.client_result_endpoint
            assert re.sub(r'\*', hostname, router.downstream_announce_endpoint) == client.upstream_announce_endpoint != self.client_announce_endpoint
            assert re.sub(r'\*', hostname, router.downstream_listen_endpoint) == client.upstream_listen_endpoint != self.client_listen_endpoint

            router.remove_router_info_file()

    def test_client_explicit_endpoint(self):
        '''Coordinated environ (w/ router): successful initialization when only client endpoints explicitly specified'''


        filename = 'server_info_filename'

        os.environ['WM_ZMQ_WRITE_INFO'] = filename

        hostname = socket.gethostname()

        with ZMQWorkManager.from_environ() as server:

            assert os.access(filename, os.F_OK)

            os.environ['WM_ZMQ_DOWNSTREAM_TASK_ENDPOINT'] = self.client_task_endpoint
            os.environ['WM_ZMQ_DOWNSTREAM_RESULT_ENDPOINT'] = self.client_result_endpoint
            os.environ['WM_ZMQ_DOWNSTREAM_ANNOUNCE_ENDPOINT'] = self.client_announce_endpoint
            os.environ['WM_ZMQ_DOWNSTREAM_LISTEN_ENDPOINT'] = self.client_listen_endpoint

            os.environ['WM_ZMQ_READ_INFO'] = filename

            router = ZMQRouter.from_environ()

            os.environ.pop('WM_ZMQ_READ_INFO', None)

            os.environ['WM_ZMQ_UPSTREAM_TASK_ENDPOINT'] = self.client_task_endpoint
            os.environ['WM_ZMQ_UPSTREAM_RESULT_ENDPOINT'] = self.client_result_endpoint
            os.environ['WM_ZMQ_UPSTREAM_ANNOUNCE_ENDPOINT'] = self.client_announce_endpoint
            os.environ['WM_ZMQ_UPSTREAM_LISTEN_ENDPOINT'] = self.client_listen_endpoint

            client = ZMQClient.from_environ()

            assert router.upstream_task_endpoint != router.downstream_task_endpoint
            assert router.upstream_result_endpoint != router.downstream_result_endpoint
            assert router.upstream_announce_endpoint != router.downstream_announce_endpoint
            assert router.upstream_listen_endpoint != router.downstream_listen_endpoint

            assert router.upstream_task_endpoint == re.sub(r'\*', hostname, server.master_task_endpoint) != self.server_task_endpoint
            assert router.upstream_result_endpoint == re.sub(r'\*', hostname, server.master_result_endpoint) != self.server_result_endpoint
            assert router.upstream_announce_endpoint == re.sub(r'\*', hostname, server.master_announce_endpoint) != self.server_announce_endpoint
            assert router.upstream_listen_endpoint == re.sub(r'\*', hostname, server.master_listen_endpoint) != self.server_listen_endpoint

            assert router.downstream_task_endpoint == client.upstream_task_endpoint == self.client_task_endpoint
            assert router.downstream_result_endpoint == client.upstream_result_endpoint == self.client_result_endpoint
            assert router.downstream_announce_endpoint == client.upstream_announce_endpoint == self.client_announce_endpoint
            assert router.downstream_listen_endpoint == client.upstream_listen_endpoint == self.client_listen_endpoint

            router.remove_router_info_file()

    def test_no_explicit_endpoints(self):
        '''Coordinated environ (w/ router): successful initialization when no endpoints explicitly specified (initialization via info files)'''

        server_filename = 'server_info_filename'
        router_filename = 'router_info_filename'

        hostname = socket.gethostname()

        os.environ['WM_ZMQ_WRITE_INFO'] = server_filename

        with ZMQWorkManager.from_environ() as server:

            os.environ['WM_ZMQ_READ_INFO'] = server_filename
            os.environ['WM_ZMQ_WRITE_INFO'] = router_filename

            router = ZMQRouter.from_environ()

            os.environ['WM_ZMQ_READ_INFO'] = router_filename

            client = ZMQClient.from_environ()

            assert router.upstream_task_endpoint == re.sub(r'\*', hostname, server.master_task_endpoint) != self.server_task_endpoint
            assert router.upstream_result_endpoint == re.sub(r'\*', hostname, server.master_result_endpoint) != self.server_result_endpoint
            assert router.upstream_announce_endpoint == re.sub(r'\*', hostname, server.master_announce_endpoint) != self.server_announce_endpoint
            assert router.upstream_listen_endpoint == re.sub(r'\*', hostname, server.master_listen_endpoint) != self.server_listen_endpoint

            assert re.sub(r'\*', hostname, router.downstream_task_endpoint) == client.upstream_task_endpoint != self.client_task_endpoint
            assert re.sub(r'\*', hostname, router.downstream_result_endpoint) == client.upstream_result_endpoint != self.client_result_endpoint
            assert re.sub(r'\*', hostname, router.downstream_announce_endpoint) == client.upstream_announce_endpoint != self.client_announce_endpoint
            assert re.sub(r'\*', hostname, router.downstream_listen_endpoint) == client.upstream_listen_endpoint != self.client_listen_endpoint

            router.remove_router_info_file()



