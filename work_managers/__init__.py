'''A system for parallel, remote execution of multiple arbitrary tasks.
Much of this, both in concept and execution, was inspired by (and in some 
cases based heavily on) the ``concurrent.futures`` package from Python 3.2,
with some simplifications and adaptations (thanks to Brian Quinlan and his
futures implementation).
'''

from core import WorkManager, WMFuture, FutureWatcher
import serial, threads, processes, zeromq
 
from serial import SerialWorkManager
from threads import ThreadsWorkManager
from processes import ProcessWorkManager
from zeromq import ZMQWorkManager, ZMQClient

_available_work_managers = {'serial': SerialWorkManager,
                            'threads': ThreadsWorkManager,
                            'processes': ProcessWorkManager,
                            'zmq': ZMQWorkManager}

import environment    
from environment import make_work_manager
