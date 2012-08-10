from nose.tools import raises, nottest, timed

class ExceptionForTest(Exception):
    pass

def will_succeed():
    return True

def will_fail():
    raise ExceptionForTest('failed as expected')

def fn_interrupted():
    raise KeyboardInterrupt

def will_hang():
    import sys,time
    time.sleep(sys.maxint)
    
def will_busyhang():
    while True:
        pass
    
def will_busyhang_uninterruptible():
    while True:
        try:
            pass
        except KeyboardInterrupt:
            pass
    
def identity(x):
    return x

def random_int(seed=None):
    import random,sys
    
    if seed is not None:
        random.seed(seed)
    
    return random.randint(0,sys.maxint)

class CommonWorkManagerTests:
    MED_TEST_SIZE=256
    
    def test_submit(self):
        future = self.work_manager.submit(will_succeed)
            
    def test_submit_many(self):
        futures = self.work_manager.submit([(will_succeed,(),{}) for i in xrange(self.MED_TEST_SIZE)])
            
    def test_as_completed(self):
        input = set(xrange(self.MED_TEST_SIZE))
        futures = [self.work_manager.submit(identity, i) for i in xrange(self.MED_TEST_SIZE)]
        output = set(future.get_result() for future in self.work_manager.as_completed(futures))
        assert input == output
        
    def test_wait_any(self):
        input = set(xrange(self.MED_TEST_SIZE))
        futures = [self.work_manager.submit(identity, i) for i in xrange(self.MED_TEST_SIZE)]
        output = self.work_manager.wait_any(futures).get_result()
        assert output in input
        
    def test_wait_all(self):
        input = set(xrange(self.MED_TEST_SIZE))
        futures = [self.work_manager.submit(identity, i) for i in xrange(self.MED_TEST_SIZE)]
        output = set(future.get_result() for future in self.work_manager.wait_all(futures))
        assert input == output
        
    @raises(ExceptionForTest)
    def test_exception_raise(self):
        future = self.work_manager.submit(will_fail)
        future.get_result()
        
    def test_exception_retrieve(self):
        future = self.work_manager.submit(will_fail)
        exc = future.get_exception()
        assert exc.args[0] == 'failed as expected'


class CommonParallelTests:
    def test_random_seq(self):
        tasks = [(random_int,(),{}) for n in xrange(self.MED_TEST_SIZE)]
        futures = self.work_manager.submit_many(tasks)
        self.work_manager.wait_all(futures)
        result_list = [future.get_result() for future in futures]
        result_set = set(result_list)
        assert len(result_list) == len(result_set)
        
    def test_random_seq_improper_seeding(self):
        tasks = [(random_int,(1979,),{}) for n in xrange(self.MED_TEST_SIZE)]
        futures = self.work_manager.submit_many(tasks)
        self.work_manager.wait_all(futures)
        result_list = [future.get_result() for future in futures]
        result_set = set(result_list)
        
        assert len(result_list) != len(result_set)
    