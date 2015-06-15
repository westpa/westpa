from core import ZMQWMError, ZMQWMTimeout, ZMQWMEnvironmentError, ZMQCore
from node import ZMQNode
from worker import ZMQWorker
from work_manager import ZMQWorkManager

import atexit
for cls in ZMQCore, ZMQNode, ZMQWorkManager, ZMQWorker:
    atexit.register(cls.remove_ipc_endpoints)