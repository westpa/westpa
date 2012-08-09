'''A system for parallel, remote execution of multiple arbitrary tasks.
Much of this, both in concept and execution, was inspired by (and in some 
cases based heavily on) the ``concurrent.futures`` package from Python 3.2,
with some simplifications and adaptations (thanks to Brian Quinlan and his
futures implementation).
'''

import cPickle as pickle
__metaclass__ = type

import logging
import sys, time, uuid, threading, signal
from collections import deque
from contextlib import contextmanager
log = logging.getLogger(__name__)

#class WEMDWorkManager:
class WorkManager:
    MODE_MASTER = 1
    MODE_WORKER = 2

    def __init__(self):
        self.mode = None
        self.prior_sigint_handler = None
        
    def __enter__(self):
        self.startup()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_traceback):
        self.shutdown()
        return False

    def sigint_handler(self, signum, frame):
        self.shutdown()
        if self.prior_sigint_handler in (signal.SIG_IGN,None):
            pass
        elif self.prior_sigint_handler == signal.SIG_DFL:
            raise KeyboardInterrupt
        else:
            self.prior_sigint_handler(signum, frame)
        
    def install_sigint_handler(self):
        self.prior_sigint_handler = signal.signal(signal.SIGINT, self.sigint_handler)

                                    
    def startup(self):
        '''Perform any necessary startup work, such as spawning clients, and return either MODE_MASTER or MODE_WORKER 
        depending on whether this work manager is a master (capable of distributing work) or a worker (capable only
        of performing work distributed by a master).'''
        self.mode = self.MODE_MASTER
        return self.MODE_MASTER
                                            
    def shutdown(self):
        '''Cleanly shut down any active workers.'''
        pass
        
    def submit(self, fn, *args, **kwargs):
        '''Submit a task to the work manager, returning a `WMFuture` object representing the pending
        result. ``fn(*args,**kwargs)`` will be executed by a worker, and the return value assigned as the
        result of the returned future.  The function ``fn`` and all arguments must be picklable; note
        particularly that off-path modules (like the system module and any active plugins) are not
        picklable unless pre-loaded in the worker process (i.e. prior to forking the master).''' 
        raise NotImplementedError
    
    def submit_many(self, tasks):
        '''Submit a set of tasks to the work manager, returning a list of `WMFuture` objects representing
        pending results. Each entry in ``tasks`` should be a triple (fn, args, kwargs), which will result in
        fn(*args, **kwargs) being executed by a worker. The function ``fn`` and all arguments must be
        picklable; note particularly that off-path modules are not picklable unless pre-loaded in the worker
        process.'''
        
        return [self.submit(fn,*args,**kwargs) for (fn,args,kwargs) in tasks]
    
    def as_completed(self, futures):
        '''Return a generator which yields results from the given ``futures`` as they become
        available.'''
        pending = set(futures)
        
        # See which futures have results, and install a watcher on those that do not
        with WMFuture.all_acquired(pending):
            completed = {future for future in futures if future.done}
            pending -= completed
            watcher = FutureWatcher(pending, threshold=1)
    
        # Yield available results immediately
        for future in completed:
            yield future
        del completed
        
        # Wait on any remaining results
        while pending:
            watcher.wait()
            completed = watcher.reset()
            for future in completed:
                yield future
                pending.remove(future)
    
    def wait_any(self, futures):
        '''Wait on any of the given ``futures`` and return the first one which has a result available.
        If more than one result is or becomes available simultaneously, any completed future may be returned.'''
        pending = set(futures)
        with WMFuture.all_acquired(pending):
            completed = {future for future in futures if future.done}
            
            if completed:
                # If any futures are complete, then we don't need to do anything else
                return completed.pop()
            else:
                # Otherwise, we need to install a watcher
                watcher = FutureWatcher(futures, threshold = 1)
        
        watcher.wait()
        completed = watcher.reset()
        return completed.pop()        
            
    def wait_all(self, futures):
        '''A convenience function which waits on all the given ``futures`` in order.  This function returns
        the same ``futures`` as submitted to the function as a list, indicating the order in which waits
        occurred.'''
        futures = list(futures)
        results = []
        for future in futures:
            results.append(future.result)
        return futures
            
class FutureWatcher:
    '''A device to wait on multiple results and/or exceptions with only one lock.'''
    
    def __init__(self, futures, threshold = 1):
        self.event = threading.Event()
        self.lock = threading.RLock()
        self.threshold = threshold
        self.completed = []
        
        for future in futures:
            future._add_watcher(self)
        
    def signal(self, future):
        '''Signal this watcher that the given future has results available. If this 
        brings the number of available futures above signal_threshold, this watcher's
        event object will be signalled as well.'''
        with self.lock:
            self.completed.append(future)
            if len(self.completed) >= self.threshold:
                self.event.set()
                
    def wait(self):
        '''Wait on one or more futures.'''
        return self.event.wait()
            
    def reset(self):
        '''Reset this watcher's list of completed futures, returning the list of completed futures
        prior to resetting it.''' 
        with self.lock:
            self.event.clear()
            completed = self.completed
            self.completed = []
            return completed
                    
class WMFuture:
    '''A "future", representing work which has been dispatched for completion asynchronously.'''
    
    @staticmethod
    @contextmanager
    def all_acquired(futures):
        '''Context manager to acquire all locks on the given ``futures``. Primarily for internal use.'''
        futures = list(futures)
        for future in futures:
            future._condition.acquire()
            
        yield # to contents of "with" block
        
        for future in futures:
            future._condition.release()
    
    def __init__(self, task_id=None):
        self.task_id = task_id or uuid.uuid4()

        self._condition = threading.Condition()
        self._done = False
        self._result = None
        self._exception = None
        self._traceback = None

        # a set of Events representing who is waiting on results from this future
        # this set will be cleared after the result is updated and watchers are notified        
        self._watchers = set()
        
        # a set of functions that will be called with this future as an argument when it is updated with a
        # result. This list will be cleared after the result is updated and all callbacks invoked
        self._update_callbacks = []  
                        
    def __repr__(self):
        return '<WMFuture 0x{id:x}: {self.task_id!s}>'.format(id=id(self), self=self)
    
    def __hash__(self):
        return hash(self.task_id)

    def _notify_watchers(self):
        '''Notify all watchers that this future has been updated, then deletes the list of update watchers.'''
        with self._condition:
            assert self._done
            for watcher in self._watchers:
                watcher.signal(self)
            self._watchers.clear()

    def _invoke_callbacks(self):
        '''Invoke all callbacks which have been registered on this future. Exceptions in callbacks
        will be logged and ignored.'''
        with self._condition:
            for callback in self._update_callbacks:
                try:
                    callback(self)
                except Exception:
                    # This may need to be a simple print to stderr, depending on the locking
                    # semantics of the logger.
                    log.exception('ignoring exception in result callback')
            del self._update_callbacks
            self._update_callbacks = []
    
    def _add_watcher(self, watcher):
        '''Add the given update watcher  to the internal list of watchers. If a result is available,
        returns immediately without updating the list of watchers.'''
        with self._condition:
            if self._done:
                watcher.signal(self)
                return
            else:
                self._watchers.add(watcher)
                
    def _add_callback(self, callback):
        '''Add the given update callback to the internal list of callbacks. If a result is available,
        invokes the callback immediately without updating the list of callbacks.'''
        with self._condition:
            if self._done:
                try:
                    callback(self)
                except Exception:
                    log.exception('ignoring exception in result callback')
            else:
                self._update_callbacks.append(callback)
                
    def _set_result(self, result):
        '''Set the result of this future to the given value, invoke on-completion callbacks, and notify
        watchers.'''
        with self._condition:
            self._result = result
            self._done = True
            self._condition.notify_all()
            self._invoke_callbacks()
            self._notify_watchers()
        
    def _set_exception(self, exception, traceback=None):
        '''Set the exception of this future to the given value, invoke on-completion callbacks, and notify
        watchers.'''
        
        with self._condition:
            self._exception = exception
            self._traceback = traceback
            self._done = True
            self._condition.notify_all()
            self._invoke_callbacks()
            self._notify_watchers()

    def get_result(self):
        '''Get the result associated with this future, blocking until it is available.'''
        with self._condition:
            if self._done:
                if self._exception:
                    if isinstance(self._traceback, basestring):
                        if self._traceback:
                            log.error('uncaught exception in remote function\n{}'.format(self._traceback))
                        raise self._exception
                    else:
                        raise self._exception, None, self._traceback                    
                else:
                    return self._result
            else:
                self._condition.wait()
                assert self._done
                if self._exception:
                    if isinstance(self._traceback, str):
                        log.error('uncaught exception in remote function\n{}'.format(self._traceback))
                        raise self._exception
                    else:
                        raise self._exception, None, self._traceback
                else:
                    return self._result
    result = property(get_result, None, None, get_result.__doc__)
    
    def wait(self):
        '''Wait until this future has a result or exception available.'''
        with self._condition:
            if self._done:
                return
            else:
                self._condition.wait()
                assert self._done
                return    

    def get_exception(self):
        '''Get the exception associated with this future, blocking until it is available.'''
        with self._condition:
            if self._done:
                return self._exception
            else:
                self._condition.wait()
                assert self._done
                return self._exception
    exception = property(get_exception, None, None, get_exception.__doc__)
    
    def get_traceback(self):
        '''Get the traceback object associated with this future, if any.'''
        with self._condition:
            if self._returned:
                return self._traceback
            else:
                self._condition.wait()
                assert self._done
                return self._traceback
    traceback = property(get_traceback, None, None, get_traceback.__doc__)
    
    def is_done(self):
        'Indicates whether this future is done executing (may block if this future is being updated).'
        with self._condition:
            return self._done
    done = property(is_done, None, None, is_done.__doc__)
    
# end class WMFuture

    
import serial, threads
import ops