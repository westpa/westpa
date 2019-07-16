# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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


import sys, logging, threading, multiprocessing
import queue
from . import WorkManager, WMFuture
import work_managers

log = logging.getLogger(__name__)

class Task:        
    def __init__(self, fn, args, kwargs, future):
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.future = future
        
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as e:
            self.future._set_exception(e, sys.exc_info()[2])
        else:
            self.future._set_result(result)
            
ShutdownSentinel = object()

class ThreadsWorkManager(WorkManager):
    '''A work manager using threads.'''
    @classmethod
    def from_environ(cls, wmenv=None):
        if wmenv is None:
            wmenv = work_managers.environment.default_env 
        return cls(wmenv.get_val('n_workers', multiprocessing.cpu_count(), int))        
    
    def __init__(self, n_workers = None):
        super(ThreadsWorkManager,self).__init__()
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self.workers = []
        self.task_queue = queue.Queue()
        
    def runtask(self, task_queue):
        while True:
            task = task_queue.get()
            if task is ShutdownSentinel:
                return
            else:
                task.run()

    def submit(self, fn, args=None, kwargs=None):
        ft = WMFuture()
        task = Task(fn, args if args is not None else (), kwargs if kwargs is not None else {}, ft)
        self.task_queue.put(task)
        return ft
                
    def startup(self):
        if not self.running:
            self.running = True
            self.workers = [threading.Thread(target=self.runtask, args=[self.task_queue], name='worker-{:d}'.format(i)) 
                            for i in range(0, self.n_workers)]
            for thread in self.workers:
                log.debug('starting thread {!r}'.format(thread))
                thread.start()
        
    def shutdown(self):
        if self.running:
            # Put one sentinel on the queue per worker, then wait for threads to terminate
            for i in range(0, self.n_workers):
                self.task_queue.put(ShutdownSentinel)
            for thread in self.workers:
                thread.join()
            self.running = False
    
    