'''
Created on Jun 11, 2015

@author: mzwier
'''

import time

from work_managers.zeromq import ZMQMaster
from . import ZMQTestBase, SETUP_WAIT, TEARDOWN_WAIT

class TestZMQMaster(ZMQTestBase):
    
    def setUp(self):
        super(TestZMQMaster,self).setUp()
        
        self.upstream_ann_endpoint = self.make_endpoint()
        self.upstream_task_endpoint = self.make_endpoint()
        self.upstream_result_endpoint = self.make_endpoint()
        
        self.downstream_ann_endpoint = self.make_endpoint()
        self.downstream_rr_endpoint = self.make_endpoint()
        
        self.test_master = ZMQMaster(self.upstream_task_endpoint, self.upstream_result_endpoint, self.upstream_ann_endpoint)
        self.test_master.downstream_ann_endpoints.append(self.downstream_ann_endpoint)
        self.test_master.downstream_rr_endpoints.append(self.downstream_rr_endpoint)
        self.test_master.startup()
        
        time.sleep(SETUP_WAIT)

    def tearDown(self):
        time.sleep(TEARDOWN_WAIT)
        
        self.test_master.signal_shutdown()
        self.test_master.comm_thread.join()
        
        super(TestZMQMaster,self).tearDown()
        

    def test_internal_shutdown(self):
        #'''master shuts down on inproc signal'''
        self.test_master.signal_shutdown()
        self.test_master.join()
        assert not self.test_master.comm_thread.is_alive()