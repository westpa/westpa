'''Routines for configuring the work manager environment'''

import os, sys, json, ConfigParser, multiprocessing

def make_work_manager(default='processes'):
    '''Using cues from the environment, instantiate a pre-configured work manager.'''
    
    work_manager_name = os.environ.get('WWMGR_WORK_MANAGER', default)
    
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
        
def get_worker_count():
    '''Get worker count from the environment'''
    return int(os.environ.get('WWMGR_N_WORKERS', multiprocessing.cpu_count()))

def get_server_info_filename():
    '''Get server info filename from the environment'''
    return os.path.abspath(os.environ.get('WWMGR_SERVER_INFO', 'wwmgr_server_info_{:d}'.format(os.getpid())))

def write_server_info(server_info, filename=None):
    json.dump(server_info, open(filename or get_server_info_filename(), 'wt'))

def read_server_info(filename=None):
    return json.load(open(filename or get_server_info_filename, 'rt'))