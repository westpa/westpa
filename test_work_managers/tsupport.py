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
