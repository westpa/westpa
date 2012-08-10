'''Routines for configuring the work manager environment'''

import os, sys, json, ConfigParser, multiprocessing

def make_work_manager(default='processes'):
    '''Using cues from the environment, instantiate a pre-configured work manager.'''
    
    work_manager_name = os.environ.get('WWMGR_WORK_MANAGER', default).lower()
    
    if work_manager_name == 'serial':
        from work_managers.serial import SerialWorkManager
        return SerialWorkManager.from_environ()
    elif work_manager_name == 'threads':
        from work_managers.threads import ThreadsWorkManager
        return ThreadsWorkManager.from_environ()
    elif work_manager_name == 'processes':
        from work_managers.processes import ProcessWorkManager
        return ProcessWorkManager.from_environ()
    elif work_manager_name in ('zeromq', 'zmq'):
        from work_managers.zeromq import ZMQWorkManager
        return ZMQWorkManager.from_environ()
    else:
        raise ValueError('unknown work manager {!r}'.format(work_manager_name))
    
def make_client():
    work_manager_name = os.environ.get('WWMGR_WORK_MANAGER', 'zmq').lower()
    if work_manager_name not in ('zmq', 'zeromq'):
        raise ValueError('work manager {} does not offer a client mode'.format(work_manager_name))
    else:
        from work_managers.zeromq import ZMQClient
        return ZMQClient.from_environ()
        
def get_worker_count():
    '''Get worker count from the environment'''
    return int(os.environ.get('WWMGR_N_WORKERS', multiprocessing.cpu_count()))

