'''
Created on May 29, 2015

@author: mzwier
'''

from __future__ import division, print_function; __metaclass__ = type

MSG_TASK_REQUEST = 'task_request'
MSG_TASK_RESULT = 'task_result'

MSG_MASTER_Q_ALIVE = 'is_master_alive'
MSG_MASTER_ACK_ALIVE = 'master_is_alive'

MSG_WORKER_Q_EXISTS = 'who_exists'
MSG_WORKER_ACK_EXISTS = 'I_exist'

MSG_SHUTDOWN = 'shutdown'
MSG_SHUTDOWN_ACK = 'shutting_down'

# Every ten seconds the master requests a status report from workers.
# This also notifies workers that the master is still alive
DEFAULT_STATUS_POLL = 10 

# If we haven't heard from the master or a worker (as appropriate) in these
# amounts of time, we assume a crash and shut down.
MASTER_CRASH_TIMEOUT = DEFAULT_STATUS_POLL * 6
WORKER_CRASH_TIMEOUT = DEFAULT_STATUS_POLL * 3

def make_message(message_code, *args, **kwargs):
    return (message_code, args, kwargs)

import gevent
import os, sys, signal
import zmq.green as zmq

signames = {val:name for name, val in reversed(sorted(signal.__dict__.items()))
            if name.startswith('SIG') and not name.startswith('SIG_')}


class ZMQCore:
    def __init__(self):
        self.context = None
        self.rr_endpoint = None
        self.ann_endpoint = None
        self.rr_socket = None
        self.ann_socket = None
        
        
    def logmessage(self, message):
        print('<{!s}> {!s}'.format(self.__class__.__name__, message))
    
    def shutdown_handler(self, signal=None):
        if signal is None:
            self.logmessage('shutting down')
        else:
            self.logmessage('shutting down on signal {!s}'.format(signames.get(signal,signal)))
        sys.exit()
    
    def install_signal_handlers(self, signals = None):
        if not signals:
            signals = {signal.SIGINT, signal.SIGQUIT, signal.SIGTERM}
        
        for sig in signals:
            gevent.signal(sig, self.shutdown_handler, sig)

    