'''Script to initialize a ZMQRouter with environment/command line args'''
from __future__ import division

import os, sys, argparse, traceback, logging

log = logging.getLogger(__name__)

import logging.config

from work_managers.zeromq import ZMQRouter
import work_managers.environment

parser = argparse.ArgumentParser('w_zmqrouter', 'startup a ZMQ Router')

group = parser.add_argument_group('general options')
egroup = group.add_mutually_exclusive_group()
egroup.add_argument('--quiet', dest='verbosity', action='store_const', const='quiet',
                     help='emit only essential information')
egroup.add_argument('--verbose', dest='verbosity', action='store_const', const='verbose',
                     help='emit extra information')
egroup.add_argument('--debug', dest='verbosity', action='store_const', const='debug',
                     help='enable extra checks and emit copious information')

ZMQRouter.add_wm_args(parser)

args = parser.parse_args()

verbosity = args.verbosity
process_name = __name__

##Configure logging
        
logging_config = {'version': 1, 'incremental': False,
                  'formatters': {'standard': {'format': '-- %(levelname)-8s [%(name)s] -- %(message)s'},
                  'debug':    {'format': '''\-- %(levelname)-8s %(asctime)24s PID %(process)-12d TID %(thread)-20d
                  from logger "%(name)s" at location %(pathname)s:%(lineno)d [%(funcName)s()] :: %(message)s'''}},
                  'handlers': {'console': {'class': 'logging.StreamHandler',
                  'stream': 'ext://sys.stdout', 'formatter': 'standard'}},
                  'loggers': {'west': {'handlers': ['console'], 'propagate': False},
                  'oldtools': {'handlers': ['console'], 'propagate': False},
                  'westtools': {'handlers': ['console'], 'propagate': False},
                  'westext': {'handlers': ['console'], 'propagate': False},
                  'work_managers': {'handlers': ['console'], 'propagate': False},
                  'multiprocessing': {'handlers': ['console'], 'propagate': False}},
                  'root': {'handlers': ['console']}}
        
logging_config['loggers'][process_name] = {'handlers': ['console'], 'propagate': False}
            
if verbosity == 'debug':
    import multiprocessing
    multiprocessing.log_to_stderr(multiprocessing.SUBDEBUG)
    logging_config['root']['level'] = 5 #'DEBUG'
    logging_config['handlers']['console']['formatter'] = 'debug'
elif verbosity == 'verbose':
    logging_config['root']['level'] = 'INFO'
else:
    logging_config['root']['level'] = 'WARNING'

logging.config.dictConfig(logging_config)
logging_config['incremental'] = True


work_managers.environment.process_wm_args(args)

router = ZMQRouter.from_environ()

router.startup()







