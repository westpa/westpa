# Copyright (C) 2013 Kim Wong and Matthew C. Zwier and Lillian T. Chong
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


"""
A work manager which uses MPI to distribute tasks and collect results.
"""


import logging, re, threading
from collections import deque
from mpi4py import MPI

import work_managers
from work_managers import WorkManager, WMFuture

log = logging.getLogger(__name__)

class Task:
    # tasks are tuples of (task_id, function, args, keyword args)
    def __init__(self, task_id, fn, args, kwargs):
        self.task_id = task_id
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        
    def __repr__(self):
        return '<Task {self.task_id}: {self.fn!r}(*{self.args!r}, **{self.kwargs!r})>'\
               .format(self=self)
    
class MPIBase:

    def __init__(self):
        # Initialize communicator and obtain standard MPI variables
        comm = MPI.COMM_WORLD

        self.comm = comm
        self.rank = comm.Get_rank()
        self.num_procs = comm.Get_size()
        self.name = MPI.Get_processor_name()

        # Define master rank
        self.master_rank = 0

        # Define message tags for task, result, and announce
        self.task_tag = 10
        self.result_tag = 20
        self.announce_tag = 30

        # create an empty message buffer
        messages = []

    def startup(self):
        raise NotImplementedError
    
    def shutdown(self):
        raise NotImplementedError

    @property
    def is_master(self):
        '''True if this is the master process for task distribution. This is necessary, e.g., for
        MPI, where all processes start identically and then must branch depending on rank.'''
        if self.rank == self.master_rank:
            return True
        else:
            return False
    
class MPIWMServer(MPIBase):
    
    def __init__(self):
        super(MPIWMServer, self).__init__()
        
        # tasks awaiting dispatch
        self.task_queue = deque()
        
        # MPI destination ranks for tasks; exclude master_rank
        self.task_dest = deque()
        for rank in range(self.num_procs):
            if rank != self.master_rank:
                self.task_dest.append(rank)

        # futures corresponding to tasks
        self.pending_futures = dict()

    def _dispatch_loop(self):
        comm = self.comm        

        while True:

            # Dispatch as many tasks as possible before checking for shutdown
            while self.task_dest:
                try:
                    task = self.task_queue.popleft()
                    task_dest = self.task_dest.popleft()
                except IndexError:
                    break
                else:
                    comm.send(task, dest = task_dest, tag = self.task_tag )

            status = MPI.Status()
            comm.Iprobe(self.master_rank, self.announce_tag, status)
            message_tag = status.Get_tag()

            # Check for announcements
            if message_tag == self.announce_tag:
                messages = comm.recv(source = self.master_rank, tag = self.announce_tag)
                if 'shutdown' in messages:
                    log.debug('exiting _dispatch_loop()')
                    return

    def _receive_loop(self):
        comm = self.comm 

        while True:

            status = MPI.Status()
            comm.Iprobe(MPI.ANY_SOURCE, MPI.ANY_TAG, status)
            message_src = status.Get_source()
            message_tag = status.Get_tag()

            # results are tuples of (task_id, {'result', 'exception'}, value)
            if message_tag == self.result_tag:
                (task_id, result_stat, result_value) = comm.recv(source = message_src, tag = message_tag)

                ft = self.pending_futures.pop(task_id)

                if result_stat == 'exception':
                    ft._set_exception(*result_value)
# Check with Matt on what else to do for an exception
                else:
                    ft._set_result(result_value)
                    self.task_dest.append(message_src)

            # Check for announcements
            elif message_tag == self.announce_tag:
                messages = comm.recv(source = message_src, tag = message_tag)
                if 'shutdown' in messages:
                    log.debug('exiting _receive_loop()')
                    return
                
    def _make_append_task(self, fn, args, kwargs):
        ft = WMFuture()
        task_id = ft.task_id
        task = Task(task_id, fn, args, kwargs)
        self.pending_futures[task_id] = ft
        self.task_queue.append(task)
        return ft
    
    def submit(self, fn, args=None, kwargs=None):
        ft = self._make_append_task(fn, args if args is not None else [], kwargs if kwargs is not None else {})
        return ft

    def startup(self):
        # start up server threads
        server_threads = []

        self._dispatch_thread = threading.Thread(target=self._dispatch_loop)
        self._dispatch_thread.start()
        server_threads.append(self._dispatch_thread)

        self._receive_thread = threading.Thread(target=self._receive_loop)
        self._receive_thread.start()
        server_threads.append(self._receive_thread)

        self.server_threads = server_threads 

class MPIClient(MPIBase):

    def __init__(self):
        super(MPIClient,self).__init__()
        
    def _create_worker(self):
        comm = self.comm

        while True:

            status = MPI.Status()
            comm.Probe(self.master_rank, MPI.ANY_TAG, status)
            message_src = self.master_rank
            message_tag = status.Get_tag()

            # Check for available task 
            if message_tag == self.task_tag:

                task = comm.recv(source = message_src, tag = message_tag)

                try:
                    result_value = task.fn(*task.args, **task.kwargs)
                except Exception as e:
                    result_object = (task.task_id, 'exception', result_value)
                else:
                    result_object = (task.task_id, 'result', result_value)

                comm.send(result_object, dest = self.master_rank, tag = self.result_tag)

            # Check for announcements
            if message_tag == self.announce_tag:
                messages = comm.recv(source = message_src, tag = message_tag)
                if 'shutdown' in messages:
                    return

    def startup(self):
        # start up client thread
        self._worker_thread = threading.Thread(target=self._create_worker)
        self._worker_thread.start()

    def run(self):
        self._worker_thread.join()

class MPIWorkManager(MPIWMServer,MPIClient,WorkManager):
    '''A work manager using MPI.'''
    @classmethod
    def from_environ(cls, wmenv=None):
        return cls()

    def __init__(self):
        WorkManager.__init__(self)
        MPIWMServer.__init__(self)
        MPIClient.__init__(self)
        
    def startup(self):
        if self.rank == self.master_rank:
            MPIWMServer.startup(self)
        else:
            MPIClient.startup(self)
            
    def shutdown(self):
        comm = self.comm
        if self.rank == self.master_rank:
            # send 'shutdown' to client threads
            for x in self.task_dest:
                comm.send('shutdown', dest = x, tag = self.announce_tag )
            # send 'shutdown' to server threads
            for thread in self.server_threads:
                comm.send('shutdown', dest = 0, tag = self.announce_tag )

        log.info( "MPIWMServer.shutdown complete" )

