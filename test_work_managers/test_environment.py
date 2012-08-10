import os
from work_managers.environment import make_work_manager
from work_managers import SerialWorkManager, ThreadsWorkManager, ProcessWorkManager, ZMQWorkManager
from tsupport import will_succeed

class TestInstantiations:
    '''Test to see that the environment system at least selects the proper work manager and sets the
    number of workers appropriately'''
    
    sanitize_vars = ('WWMGR_WORK_MANAGER', 'WWMGR_N_WORKERS')
    
    def setUp(self):
        for varname in self.sanitize_vars:
            assert varname not in os.environ
    
    def tearDown(self):
        for varname in self.sanitize_vars:
            os.environ.pop(varname, None)
        
    def testSerial(self):
        os.environ['WWMGR_WORK_MANAGER'] = 'serial'
        work_manager = make_work_manager()
        assert isinstance(work_manager, SerialWorkManager)
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        
    def testThreads(self):
        os.environ['WWMGR_WORK_MANAGER'] = 'threads'
        os.environ['WWMGR_N_WORKERS'] = str(3)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ThreadsWorkManager)
        assert work_manager.n_workers == 3
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        
    def testProcesses(self):
        os.environ['WWMGR_WORK_MANAGER'] = 'processes'
        os.environ['WWMGR_N_WORKERS'] = str(3)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ProcessWorkManager)
        assert work_manager.n_workers == 3
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        
    def testZeroMQ(self):
        os.environ['WWMGR_WORK_MANAGER'] = 'zmq'
        os.environ['WWMGR_N_WORKERS'] = str(3)
        work_manager = make_work_manager()
        assert isinstance(work_manager, ZMQWorkManager)
        assert work_manager.internal_client.n_workers == 3
        with work_manager:
            future = work_manager.submit(will_succeed)
            future.get_result()
        