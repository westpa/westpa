# Copyright (C) 2013 Matthew C. Zwier, Joshua L. Adelman and Lillian T. Chong (see note below).
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
#
# This implementation is derived from the ``concurrent.futures``
# module of Python 3.2, by Brian Quinlan, (C) 2011 the Python Software
# Foundation. See http://docs.python.org/3/license.html for more information.

import logging
import uuid, threading, signal
from itertools import islice
from contextlib import contextmanager
log = logging.getLogger(__name__)

class WorkManager:
    '''Base class for all work managers. At a minimum, work managers must provide a 
    ``submit()`` function and a ``n_workers`` attribute (which may be a property),
    though most will also override ``startup()`` and ``shutdown()``.'''
    
    @classmethod
    def from_environ(cls, wmenv=None):
        raise NotImplementedError
    
    @classmethod
    def add_wm_args(cls, parser, wmenv=None):
        return
    
    def __repr__(self):
        return '<{classname} at 0x{id:x}>'.format(classname=self.__class__.__name__, id=id(self))
    
    def __init__(self):
        self._sigint_handler_installed = False
        self.prior_sigint_handler = None
        self.running = False
        
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
        if not self._sigint_handler_installed:
            self._sigint_handler_installed = True
            self.prior_sigint_handler = signal.signal(signal.SIGINT, self.sigint_handler)
                                    
    def startup(self):
        '''Perform any necessary startup work, such as spawning clients.'''
        self.running = True
                                            
    def shutdown(self):
        '''Cleanly shut down any active workers.'''
        self.running = False
        
    
    def run(self):
        '''Run the worker loop (in clients only).'''
        pass
        
    def submit(self, fn, args=None, kwargs=None):
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
        
        return [self.submit(fn,args,kwargs) for (fn,args,kwargs) in tasks]
    
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

    def submit_as_completed(self, task_generator, queue_size=None):
        '''Return a generator which yields results from a set of ``futures`` as they become
        available. Futures are generated by the ``task_generator``, which must return a triple of the form
        expected by ``submit``. The method also accepts an int ``queue_size`` that dictates the
        maximum number of Futures that should be pending at any given time. The default value of
        ``None`` submits all of the tasks at once.'''

        futures = [self.submit(fn,args,kwargs) for (fn,args,kwargs) in islice(task_generator, queue_size)]
        pending = set(futures)

        with WMFuture.all_acquired(pending):
            watcher = FutureWatcher(pending, threshold=1)

        while pending:
            watcher.wait()
            completed = watcher.reset()
            new_futures = [self.submit(fn,args,kwargs) for (fn,args,kwargs) in islice(task_generator, len(completed))]
            pending.update(new_futures)

            with WMFuture.all_acquired(new_futures):
                watcher.add(new_futures)

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

    @property
    def is_master(self):
        '''True if this is the master process for task distribution. This is necessary, e.g., for
        MPI, where all processes start identically and then must branch depending on rank.'''
        return True
            

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

    def add(self, futures):
        '''Add watchers to all futures in the iterable of futures.''' 
        for future in futures:
            future._add_watcher(self)


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

    def get_result(self,discard=True):
        '''Get the result associated with this future, blocking until it is available.
        If ``discard`` is true, then removes the reference to the result contained
        in this instance, so that a collection of futures need not turn into a cache of
        all associated results.'''
        with self._condition:
            if self._done:
                if self._exception:
                    if isinstance(self._traceback, str):
                        if self._traceback:
                            log.error('uncaught exception in remote function\n{}'.format(self._traceback))
                        raise self._exception
                    else:
                        raise self._exception.with_traceback(self._traceback)
            else:
                self._condition.wait()
                assert self._done
                if self._exception:
                    if isinstance(self._traceback, str):
                        log.error('uncaught exception in remote function\n{}'.format(self._traceback))
                        raise self._exception
                    else:
                        raise self._exception.with_traceback(self._traceback)
                
            result = self._result
            if discard:
                del self._result
            return result
        
    @property
    def result(self):
        return self.get_result(discard=False)
        
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
