'''
Created on May 29, 2015

@author: mzwier
'''

from __future__ import division, print_function; __metaclass__ = type


# Every ten seconds the master requests a status report from workers.
# This also notifies workers that the master is still alive
DEFAULT_STATUS_POLL = 10 

# If we haven't heard from the master or a worker (as appropriate) in these
# amounts of time, we assume a crash and shut down.
MASTER_CRASH_TIMEOUT = DEFAULT_STATUS_POLL * 6
WORKER_CRASH_TIMEOUT = DEFAULT_STATUS_POLL * 3

import logging
log = logging.getLogger(__name__)

import gevent
import sys, uuid

import signal
signames = {val:name for name, val in reversed(sorted(signal.__dict__.items()))
            if name.startswith('SIG') and not name.startswith('SIG_')}

class ZMQWMError(RuntimeError):
    '''Base class for errors related to the ZeroMQ work manager itself'''
    pass

class ZMQWMEnvironmentError(ZMQWMError):
    '''Class representing an error in the environment in which the ZeroMQ work manager is running.
    This includes such things as master/worker ID mismatches.'''
    
class ZMQWMTimeout(ZMQWMEnvironmentError):
    '''A timeout of a sort that indicatess that a master or worker has failed or never started.'''

class Message:
    ACK = 'ok'
    IDENTIFY_MASTER = 'identify_master'
    TASK_REQUEST = 'task_request'
    TASK_RESULT = 'task_result'
    SHUTDOWN = 'shutdown'
    MASTER_BEACON = 'master alive'
    
    def __init__(self, message=None, payload=None, master_id=None, worker_id=None):
        self.master_id = master_id
        self.worker_id = worker_id
        self.message = message
        self.payload = payload
        
    def __repr__(self):
        return ('<{!s} master_id={master_id!s} worker_id={worker_id!s} message={message!r} payload={payload!r}>'
                .format(self.__class__.__name__, **self.__dict__))   

class ZMQCore:
    def __init__(self):
        self.context = None
        self.rr_endpoint = None
        self.ann_endpoint = None
        self.rr_socket = None
        self.ann_socket = None
        
        self.id = uuid.uuid4()
        self.master_id = None
        
    def __repr__(self):
        return '<{!s} {!s}>'.format(self.__class__.__name__, self.id)

    def validate_message(self, message):
        '''Validate incoming message. Raises an exception if the message is improperly formatted (TypeError)
        or does not correspond to the appropriate master (ZMQWMEnvironmentError).'''
        try:
            super_validator = super(ZMQCore,self).validate_message
        except AttributeError:
            pass
        else:
            super_validator(message)
        
    def recv_message(self, socket, validate=True):
        message = socket.recv_pyobj()
        if not isinstance(message, Message):    
            raise TypeError('message is not an instance of core.Message')        
        if validate:
            self.validate_message(message)
        return message
        
    def dumpinfo(self, message):
        print('{!r} {!s}'.format(self, message))
    
    def shutdown_handler(self, signal=None):
        if signal is None:
            self.dumpinfo('shutting down')
        else:
            self.dumpinfo('shutting down on signal {!s}'.format(signames.get(signal,signal)))
        sys.exit()
    
    def install_signal_handlers(self, signals = None):
        if not signals:
            signals = {signal.SIGINT, signal.SIGQUIT, signal.SIGTERM}
        
        for sig in signals:
            gevent.signal(sig, self.shutdown_handler, sig)

    