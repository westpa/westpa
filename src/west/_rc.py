# WEST run control routines

from __future__ import division, print_function; __metaclass__ = type

import logging
log = logging.getLogger('west.rc')

import os, sys, errno
import west
from west.util.config_dict import ConfigDict
from west.util import extloader


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

class _WESTRC:
    '''A class, an instance of which is accessible as ``west.rc``, to handle global issues for WEST code,
    such as loading modules and plugins, writing output based on verbosity level, adding default command line options,
    and so on.'''
    
    # Runtime config file management
    ENV_RUNTIME_CONFIG  = 'WESTRC'
    RC_DEFAULT_FILENAME = 'west.cfg'
    
    DEFAULT_WORK_MANAGER = 'threads'
        
    def __init__(self):        
        self.verbosity = None
        self.rcfile = os.environ.get(self.ENV_RUNTIME_CONFIG) or self.RC_DEFAULT_FILENAME

        self.config = ConfigDict()
        self.process_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    
        self._system = None
        self._data_manager = None
        self._sim_manager = None
        self._we_driver = None
        self._propagator = None
        
        self.status_stream = sys.stdout
        
    def add_args(self, parser):
        group = parser.add_argument_group('general options')
        group.add_argument('-r', '--rcfile', metavar='RCFILE', dest='rcfile',
                            default=(os.environ.get(self.ENV_RUNTIME_CONFIG) or self.RC_DEFAULT_FILENAME),
                            help='use RCFILE as the WEST run-time configuration file (default: %(default)s)')
        
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--quiet', dest='verbosity', action='store_const', const='quiet',
                             help='emit only essential information')
        egroup.add_argument('--verbose', dest='verbosity', action='store_const', const='verbose',
                             help='emit extra information')
        egroup.add_argument('--debug', dest='verbosity', action='store_const', const='debug',
                            help='enable extra checks and emit copious information')
        
        group.add_argument('--version', action='version', version='WEST version %s' % west.version)
        
    @property
    def verbose_mode(self):
        return (self.verbosity in ('verbose', 'debug'))
    
    @property
    def debug_mode(self):
        return (self.verbosity == 'debug')
    
    @property
    def quiet_mode(self):
        return (self.verbosity == 'quiet')
                            
    def process_args(self, args, config_required = True):
        self.cmdline_args = args
        self.verbosity = args.verbosity
        
        if args.rcfile:
            self.rcfile = args.rcfile
        
        try:
            self.read_config()
        except IOError as e:
            if e.errno == errno.ENOENT and not config_required:
                pass
            else:
                raise
        self.config_logging()
        self.config.update_from_object(args)
                            
    def read_config(self, filename = None):
        if filename:
            self.rcfile = filename

        if 'WEST_SIM_ROOT' not in os.environ:
            sys.stderr.write('-- WARNING -- setting $WEST_SIM_ROOT to current directory ({})\n'.format(os.getcwd()))
            os.environ['WEST_SIM_ROOT'] = os.getcwd()
                                    
        self.config.read_config_file(self.rcfile) 
                    
    def config_logging(self):
        import logging.config
        logging_config = {'version': 1, 'incremental': False,
                          'formatters': {'standard': {'format': '-- %(levelname)-8s [%(name)s] -- %(message)s'},
                                         'debug':    {'format': '''\
-- %(levelname)-8s %(asctime)24s PID %(process)-12d TID %(thread)-20d
   from logger "%(name)s" 
   at location %(pathname)s:%(lineno)d [%(funcName)s()] 
   ::
   %(message)s
'''}},
                          'handlers': {'console': {'class': 'logging.StreamHandler',
                                                   'stream': 'ext://sys.stdout',
                                                   'formatter': 'standard'}},
                          'loggers': {'west': {'handlers': ['console'], 'propagate': False},
                                      'oldtools': {'handlers': ['console'], 'propagate': False},
                                      'westext': {'handlers': ['console'], 'propagate': False},
                                      'work_managers': {'handlers': ['console'], 'propagate': False},
                                      'multiprocessing': {'handlers': ['console'], 'propagate': False}},
                          'root': {'handlers': ['console']}}
        
        logging_config['loggers'][self.process_name] = {'handlers': ['console'], 'propagate': False}
            
        if self.verbosity == 'debug':
            import multiprocessing
            multiprocessing.log_to_stderr(multiprocessing.SUBDEBUG)
            logging_config['root']['level'] = 5 #'DEBUG'
            logging_config['loggers']['multiprocessing']['level'] = 'DEBUG'
            logging_config['handlers']['console']['formatter'] = 'debug'
        elif self.verbosity == 'verbose':
            logging_config['root']['level'] = 'INFO'
        else:
            logging_config['root']['level'] = 'WARNING'

        logging.config.dictConfig(logging_config)
        logging_config['incremental'] = True
        
    def pstatus(self, *args, **kwargs):
        fileobj = kwargs.pop('file', self.status_stream)
        if kwargs.get('termonly', False) and not fileobj.isatty():
            return
        if self.verbosity != 'quiet':
            print(*args, file=fileobj, **kwargs)
            
    def pstatus_term(self, *args, **kwargs):
        fileobj = kwargs.pop('file', self.status_stream)
        if fileobj.isatty() and self.verbosity != 'quiet':
            print(*args, file=fileobj, **kwargs)
        
    def pflush(self):
        for stream in (self.status_stream, sys.stdout, sys.stderr):
            try:
                stream.flush()
            except AttributeError:
                pass
            
    def new_sim_manager(self, work_manager=None):
        drivername = self.config.get('drivers.sim_manager', 'default')
        if drivername.lower() == 'default':
            from west.sim_manager import WESimManager
            sim_manager = WESimManager(work_manager)
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path')
            sim_manager = extloader.get_object(drivername,pathinfo)(work_manager)
        log.debug('loaded simulation manager {!r}'.format(sim_manager))
        return sim_manager
        
    def get_sim_manager(self, work_manager=None):
        if self._sim_manager is None:
            self._sim_manager = self.new_sim_manager(work_manager)
        return self._sim_manager

    def new_data_manager(self):
        drivername = self.config.get('drivers.data_manager', 'hdf5')
        if drivername.lower() in ('hdf5', 'default'):
            data_manager = west.data_manager.WESTDataManager()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path', default=None)
            data_manager = extloader.get_object(drivername, pathinfo)()
        log.debug('loaded data manager: {!r}'.format(data_manager))
        return data_manager
        
    def get_data_manager(self):
        if self._data_manager is None:
            self._data_manager = self.new_data_manager()
        return self._data_manager
    
    def new_we_driver(self):
        drivername = self.config.get('drivers.we_driver', 'default')
        if drivername.lower() == 'default':
            we_driver = west.we_driver.WEDriver()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path', default=None)
            we_driver = extloader.get_object(drivername, pathinfo)()
        log.debug('loaded WE algorithm driver: {!r}'.format(we_driver))
        return we_driver
        
    def get_we_driver(self):
        if self._we_driver is None:
            self._we_driver = self.new_we_driver()
        return self._we_driver    
    
    def new_propagator(self):
        drivername = self.config.require('drivers.propagator')
        if drivername.lower() == 'executable':
            import west.propagators.executable
            propagator = west.propagators.executable.ExecutablePropagator()
        else:
            pathinfo = self.config.get_pathlist('drivers.module_path', default=None)
            propagator = extloader.get_object(drivername, pathinfo)()
        log.debug('loaded propagator {!r}'.format(propagator))
        return propagator
        
    def get_propagator(self):
        if self._propagator is None:
            self._propagator = self.new_propagator()
        return self._propagator
            
    def new_system_driver(self):
        sysdrivername = self.config.require('system.system_driver')
        log.info('loading system driver %r' % sysdrivername)
        pathinfo = self.config.get_pathlist('system.module_path', default=None)
        pathinfo.append(os.environ.get('WEST_SIM_ROOT', '.'))

        system = extloader.get_object(sysdrivername, pathinfo)()
        system.initialize()
        log.debug('loaded system driver {!r}'.format(system))        
        return system
    
    def get_system_driver(self):
        if self._system is None:
            self._system = self.new_system_driver()
        return self._system
    
    propagator = property(get_propagator)
    we_driver = property(get_we_driver)
    system = property(get_system_driver)
    data_manager = property(get_data_manager)
    sim_manager = property(get_sim_manager)
