'''A system for parallel, remote execution of multiple arbitrary tasks.
Much of this, both in concept and execution, was inspired by (and in some
cases based heavily on) the ``concurrent.futures`` package from Python 3.2,
with some simplifications and adaptations (thanks to Brian Quinlan and his
futures implementation).
'''

import logging

from .core import WorkManager, WMFuture, FutureWatcher  # noqa


# Import core work managers, which should run most everywhere that
# Python does
from . import serial, threads, processes  # noqa
from .serial import SerialWorkManager
from .threads import ThreadsWorkManager
from .processes import ProcessWorkManager

log = logging.getLogger(__name__)

_available_work_managers = {'serial': SerialWorkManager, 'threads': ThreadsWorkManager, 'processes': ProcessWorkManager}

# Import ZeroMQ work manager if available
try:
    from . import zeromq  # noqa
    from .zeromq import ZMQWorkManager
except ImportError:
    log.info('ZeroMQ work manager not available')
    log.debug('traceback follows', exc_info=True)
else:
    _available_work_managers['zmq'] = ZMQWorkManager

# Import MPI work manager if available
try:
    from . import mpi  # noqa
    from .mpi import MPIWorkManager
except ImportError:
    log.info('MPI work manager not available')
    log.debug('traceback follows', exc_info=True)
else:
    _available_work_managers['mpi'] = MPIWorkManager

from . import environment  # noqa
from .environment import make_work_manager  # noqa


__all__ = [
    'serial',
    'threads',
    'processes',
    'SerialWorkManager',
    'ThreadsWorkManager',
    'ProcessWorkManager',
    'environment',
    'make_work_manager',
]
