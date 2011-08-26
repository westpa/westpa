# WEMD run control routines

from __future__ import division, print_function; __metaclass__ = type

import logging
log = logging.getLogger('wemd.rc')

import os, sys, argparse
import wemd
from wemd.util.config_dict import ConfigDict
from wemd.util import extloader


def lazy_loaded(backing_name, loader, docstring = None):
    def getter(self):
        obj = getattr(self, backing_name, None)
        if obj is None:
            obj = loader()
            setattr(self,backing_name,obj)
        return obj
    def setter(self, val):
        setattr(self,backing_name,val)
    def deleter(self):
        delattr(self,backing_name)
        setattr(self,backing_name,None)
    return property(getter, setter, deleter, docstring)

class _WEMDRC:
    # Runtime config file management
    ENV_RUNTIME_CONFIG  = 'WEMDRC'
    RC_DEFAULT_FILENAME = 'wemd.cfg'
    
    # Exit codes
    EX_SUCCESS           = 0
    EX_ERROR             = 1
    EX_USAGE_ERROR       = 2
    EX_ENVIRONMENT_ERROR = 3
    EX_RTC_ERROR         = 4
    EX_DB_ERROR          = 5
    EX_STATE_ERROR       = 6
    EX_CALC_ERROR        = 7
    EX_COMM_ERROR        = 8
    EX_DATA_ERROR        = 9
    EX_EXCEPTION_ERROR   = 10
    
    _ex_names = dict((code, name) for (name, code) in locals().iteritems()
                     if name.startswith('EX_'))
    
    @classmethod
    def get_exit_code_name(cls, code):
        return cls._ex_names.get(code, 'error %d' % code)

    @classmethod
    def common_arg_parser(cls,prog=None,description=None):
        '''Return an argparse.ArgumentParser pre-loaded with --rcfile, --quiet, --verbose, --debug,
        and --version'''
        parser = argparse.ArgumentParser(prog=prog,description=description)
        parser.add_argument('-r', '--rcfile', metavar='RCFILE', dest='rcfile',
                            default=(os.environ.get(cls.ENV_RUNTIME_CONFIG) or cls.RC_DEFAULT_FILENAME),
                            help='use RCFILE as the WEMD run-time configuration file (default: %(default)s)')
        parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                            help='emit only essential information')
        parser.add_argument('--verbose', dest='verbose_mode', action='store_true',
                            help='emit extra information')
        parser.add_argument('--debug', dest='debug_mode', action='store_true',
                            help='enable extra checks and emit copious information')
        parser.add_argument('--version', action='version', version='WEMD version %s' % wemd.version)
        return parser
    
    def __init__(self):        
        self.quiet_mode = None
        self.verbose_mode = None
        self.debug_mode = None
        self.rcfile = os.environ.get(self.ENV_RUNTIME_CONFIG) or self.RC_DEFAULT_FILENAME

        self.config = ConfigDict()
        self.process_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
                            
    def process_common_args(self, args):
        self.cmdline_args = args
        
        self.quiet_mode = args.quiet_mode
        self.verbose_mode = args.verbose_mode
        self.debug_mode = args.debug_mode
        if args.rcfile:
            self.rcfile = args.rcfile
        
        self.read_config()
        self.config_logging()
        self.config.update_from_object(args)
                    
    def read_config(self, filename = None):
        if filename:
            self.rcfile = filename            
        self.config.read_config_file(self.rcfile)
                    
    def config_logging(self):
        import logging.config
        logging_config = {'version': 1, 'incremental': False,
                          'formatters': {'standard': {'format': '  -- %(levelname)-8s -- %(message)s'},
                                         'debug':    {'format': '''\
          -- %(levelname)-8s %(asctime)24s PID %(process)-12d TID %(thread)-20d 
             %(pathname)s:%(lineno)d [%(funcName)s()] 
               %(message)s'''}},
                          'handlers': {'console': {'class': 'logging.StreamHandler',
                                                   'stream': 'ext://sys.stdout',
                                                   'formatter': 'standard'}},
                          'loggers': {'wemd': {'handlers': ['console'], 'propagate': False},
                                      'wemdtools': {'handlers': ['console'], 'propagate': False},
                                      'wemd_cli': {'handlers': ['console'], 'propagate': False}},
                          'root': {'handlers': ['console']}}
        
        logging_config['loggers'][self.process_name] = {'handlers': ['console'], 'propagate': False}
            
        if self.debug_mode:
            logging_config['root']['level'] = 'DEBUG'
            logging_config['handlers']['console']['formatter'] = 'debug'
        elif self.verbose_mode:
            logging_config['root']['level'] = 'INFO'
        else:
            logging_config['root']['level'] = 'WARNING'

        logging.config.dictConfig(logging_config)
        logging_config['incremental'] = True
        
    def pstatus(self, *args, **kwargs):
        if self.quiet_mode:
            print(*args, **kwargs)
        
    def get_sim_manager(self):
        drivername = self.config.get('drivers.sim_manager', 'default')
        if drivername.lower() == 'default':
            from wemd.sim_manager import WESimManager
            return WESimManager()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path')
            return extloader.get_object(drivername,pathinfo)()
        
    def get_data_manager(self):
        drivername = self.config.get('drivers.data_manager', 'hdf5')
        if drivername.lower() in ('hdf5', 'default'):
            data_manager = wemd.data_manager.WEMDDataManager()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path', default=None)
            data_manager = extloader.get_object(drivername, pathinfo)()
        log.debug('loaded data manager: {!r}'.format(data_manager))
        return data_manager

    def get_we_driver(self):
        drivername = self.config.get('drivers.we_driver', 'default')
        if drivername.lower() == 'default':
            we_driver = wemd.we_driver.WEMDWEDriver()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path', default=None)
            we_driver = extloader.get_object(drivername, pathinfo)()
        log.debug('loaded WE algorithm driver: {!r}'.format(we_driver))
        return we_driver

    def get_work_manager(self):
        drivername = self.config.get('args.work_manager_name')
        if not drivername:
            drivername = self.config.get('drivers.work_manager', 'threads')
        if drivername.lower() == 'serial':
            import wemd.work_managers.serial
            work_manager = wemd.work_managers.serial.SerialWorkManager()
        elif drivername.lower() == 'processes':
            import wemd.work_managers.processes
            work_manager = wemd.work_managers.processes.ProcessWorkManager()                    
        elif drivername.lower() in ('threads', 'default'):
            import wemd.work_managers.threads
            work_manager = wemd.work_managers.threads.ThreadedWorkManager()
        elif drivername.lower() in ('zmq', 'zeromq'):
            import wemd.work_managers.zeromq
            work_manager = wemd.work_managers.zeromq.ZMQWorkManager()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path', default=None)
            work_manager = extloader.get_object(drivername, pathinfo)()
        log.debug('loaded work manager: {!r}'.format(work_manager))
        return work_manager
    
    def get_propagator(self):
        drivername = self.config.require('drivers.propagator')
        if drivername.lower() == 'executable':
            import wemd.propagators.executable
            propagator = wemd.propagators.executable.ExecutablePropagator()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path', default=None)
            propagator = extloader.get_object(drivername, pathinfo)()
        log.debug('loaded propagator {!r}'.format(propagator))
        return propagator
    
    def get_system_driver(self):
        sysdrivername = self.config.require('system.system_driver')
        log.info('loading system driver %r' % sysdrivername)
        pathinfo = self.config.get_pathlist('system.module_path', default=None)        
        system = extloader.get_object(sysdrivername, pathinfo)()
        log.debug('loaded system driver {!r}'.format(system))
        
        log.debug('initializing system driver')
        system.initialize()

        return system
