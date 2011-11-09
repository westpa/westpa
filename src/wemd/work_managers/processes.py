from __future__ import division, print_function; __metaclass__ = type

import sys, os, logging, socket, multiprocessing, threading, time, numpy
import argparse
import collections
from functools import partial
import wemd
#from wemd.work_managers import WEMDWorkManager, Task, Result
from wemd.work_managers import WEMDWorkManager, WMFuture

log = logging.getLogger(__name__)

class ProcessWorkManager(WEMDWorkManager):
    '''A work manager using the ``multiprocessing`` module.'''
    
    def __init__(self):
        super(ProcessWorkManager,self).__init__()
        self.n_workers = None
        self.pool = None
                
        
    def _result_callback(self, task_id, asyncresult):
        future = self.futures[task_id]
        with future._condition:
            try:
                result = asyncresult.get()
            except Exception as e:
                future._set_exception(e)
            else:
                future._set_result(result)
        self._remove_future(future)
        
    def parse_aux_args(self, aux_args, do_help = False):
        parser = argparse.ArgumentParser(usage='%(prog)s [NON_WORK_MANAGER_OPTIONS] [OPTIONS]',
                                         add_help=False)
        
        runtime_config = wemd.rc.config
        group = parser.add_argument_group('processes work manager options')
                
        group.add_argument('-n', '-np', '-nw', type=int, dest='n_workers', default=multiprocessing.cpu_count(),
                            help='Number of worker processes to run. (Default: %(default)s)')
        if do_help:
            parser.print_help()
            sys.exit(0)
        args, extra_args = parser.parse_known_args(aux_args)
        
        self.n_workers = runtime_config['processes_work_manager.n_workers'] = args.n_workers
        
        return extra_args
        
    def startup(self):
        self.pool = multiprocessing.Pool(self.n_workers)
        self.dispatch_queue = collections.deque
        self.result_queue = collections.deque
        
    def shutdown(self, exit_code = 0):
        self.pool.terminate()
        self.pool.join()
        self.pool = None
        self.dispatch_queue = None
        self.result_queue = None
        
    def submit(self, fn, *args, **kwargs):
        ft = WMFuture()
        self.pool.apply_async(fn, args=args, kwds=kwargs, callback=partial(self._result_callback, ft.task_id))
        self._register_future(ft)
        return ft
    
    