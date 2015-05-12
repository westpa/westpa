# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

"""WEST run control and configuration routines"""

from __future__ import division, print_function; __metaclass__ = type

import logging
log = logging.getLogger('westpa.rc')

import os, sys, errno
import westpa
from yamlcfg import YAMLConfig
from yamlcfg import YAMLSystem
from . import extloader
from work_managers import SerialWorkManager

# For making it's own mapper 
from westpa.binning import RectilinearBinMapper
    
def bins_from_yaml_dict(bin_dict):
    typename = bin_dict.pop('type')
    kwargs = bin_dict
    
    try:
        mapper_type = getattr(sys.modules['westpa.binning'], typename)
    except AttributeError:
        raise KeyError('unknown bin mapper type {!r} in config file {!r}'.format(typename))
    
    if typename == 'RectilinearBinMapper':
        boundary_lists = kwargs.pop('boundaries')
        for ilist, boundaries in enumerate(boundary_lists):
            boundary_lists[ilist] = map((lambda x: 
                                           float('inf') 
                                           if (x if isinstance(x, basestring) else '').lower() == 'inf' 
                                           else x), boundaries)
        return mapper_type(boundary_lists)
    else:
        try:
            return mapper_type(**kwargs)
        except Exception:
            log.exception('exception instantiating mapper')
            raise

def parsePCV(pc_str):
    namespace = {'math': math,
                 'numpy': numpy,
                 'inf': float('inf')}
    
    arr = numpy.array(eval(pc_str,namespace))
    if arr.ndim == 0:
        arr.shape = (1,1)
    elif arr.ndim == 1:
        arr.shape = (1,) + arr.shape 
    else:
        raise ValueError('too many dimensions')
    return arr


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

class WESTRC:
    '''A class, an instance of which is accessible as ``westpa.rc``, to handle global issues for WEST-PA code,
    such as loading modules and plugins, writing output based on verbosity level, adding default command line options,
    and so on.'''
    
    # Runtime config file management
    ENV_RUNTIME_CONFIG  = 'WESTRC'
    RC_DEFAULT_FILENAME = 'west.cfg'
        
    def __init__(self):        
        self.verbosity = None
        self.rcfile = os.environ.get(self.ENV_RUNTIME_CONFIG) or self.RC_DEFAULT_FILENAME

        self.config = YAMLConfig()
        self.process_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    
        # Crucial simulation and analysis drivers
        self._system = None
        self._data_manager = None
        self._sim_manager = None
        self._we_driver = None
        self._propagator = None
        
        self.work_manager = SerialWorkManager()
        
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
        
        group.add_argument('--version', action='version', version='WEST version %s' % westpa.version)
        
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
        self.config['args'] = {k:v for k,v in args.__dict__.iteritems() if not k.startswith('_')}
        self.process_config()
    
    def process_config(self):
        log.debug('config: {!r}'.format(self.config)) 
        sys.path.extend(self.config.get_pathlist(['west', 'drivers', 'module_path'], []))
        try:
            sys.path.append(os.environ['WEST_SIM_ROOT'])
        except KeyError:
            pass
                            
    def read_config(self, filename = None):
        if filename:
            self.rcfile = filename

        if 'WEST_SIM_ROOT' not in os.environ:
            #sys.stderr.write('-- WARNING -- setting $WEST_SIM_ROOT to current directory ({})\n'.format(os.getcwd()))
            os.environ['WEST_SIM_ROOT'] = os.getcwd()
                                    
        self.config.update_from_file(self.rcfile) 
                    
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
                                      'westpa': {'handlers': ['console'], 'propagate': False},
                                      'oldtools': {'handlers': ['console'], 'propagate': False},
                                      'westtools': {'handlers': ['console'], 'propagate': False},
                                      'westext': {'handlers': ['console'], 'propagate': False},
                                      'work_managers': {'handlers': ['console'], 'propagate': False},
                                      'py.warnings': {'handlers': ['console'], 'propagate': False}},
                          'root': {'handlers': ['console']}}
        
        logging_config['loggers'][self.process_name] = {'handlers': ['console'], 'propagate': False}
            
        if self.verbosity == 'debug':
            logging_config['root']['level'] = 5 #'DEBUG'
            logging_config['handlers']['console']['formatter'] = 'debug'
        elif self.verbosity == 'verbose':
            logging_config['root']['level'] = 'INFO'
        else:
            logging_config['root']['level'] = 'WARNING'

        logging.config.dictConfig(logging_config)
        logging_config['incremental'] = True
        logging.captureWarnings(True)
        
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
            
    def new_sim_manager(self):
        drivername = self.config.get(['west', 'drivers', 'sim_manager'], 'default')
        if drivername.lower() == 'default':
            from west.sim_manager import WESimManager
            sim_manager = WESimManager(rc=self)
        else:
            sim_manager = extloader.get_object(drivername)(rc=self)
        log.debug('loaded simulation manager {!r}'.format(sim_manager))
        return sim_manager
        
    def get_sim_manager(self):
        if self._sim_manager is None:
            self._sim_manager = self.new_sim_manager()
        return self._sim_manager

    def new_data_manager(self):
        import west
        drivername = self.config.get(['west', 'drivers', 'data_manager'], 'hdf5')
        if drivername.lower() in ('hdf5', 'default'):
            data_manager = west.data_manager.WESTDataManager()
        else:
            data_manager = extloader.get_object(drivername)(rc=self)
        log.debug('loaded data manager: {!r}'.format(data_manager))
        return data_manager
        
    def get_data_manager(self):
        if self._data_manager is None:
            self._data_manager = self.new_data_manager()
        return self._data_manager
    
    def new_we_driver(self):
        import west
        drivername = self.config.get(['west', 'drivers', 'we_driver'], 'default')
        if drivername.lower() == 'default':
            we_driver = west.we_driver.WEDriver()
        else:
            we_driver = extloader.get_object(drivername)(rc=self)
        log.debug('loaded WE algorithm driver: {!r}'.format(we_driver))
        return we_driver
        
    def get_we_driver(self):
        if self._we_driver is None:
            self._we_driver = self.new_we_driver()
        return self._we_driver    
    
    def new_propagator(self):
        drivername = self.config.require(['west', 'propagation', 'propagator'])
        if drivername.lower() == 'executable':
            import west.propagators.executable
            propagator = west.propagators.executable.ExecutablePropagator()
        else:
            propagator = extloader.get_object(drivername)(rc=self)
        log.debug('loaded propagator {!r}'.format(propagator))
        return propagator
        
    def get_propagator(self):
        if self._propagator is None:
            self._propagator = self.new_propagator()
        return self._propagator
            
    def new_system_driver(self):
        ''' 
        1) Build the state from driver > update from YAML (if exists)
        2) If the driver doesn't exist, build directly from yaml
        '''

        # This horribly written code is supposed to be a prototype for 
        # the YAML System builder

        # First step is to see if we have a driver specified
        system = None
        try: 
            sysdrivername = self.config.get(['west', 'system', 'driver']) 
            log.info('loading system driver %r' % sysdrivername)
            system_driver = extloader.get_object(sysdrivername)(rc=self)
            system_driver.initialize()
            log.debug('loaded system driver {!r}'.format(system_driver))        
            system = system_driver
        except KeyError:
            log.info("System driver not specified")
        # Second let's see if we have info in the YAML file 
        try: 
            yamloptions = self.config['west']['system']['system_options']
            log.info("Loading system options from configuration file")
            if system:
                system = self.update_from_yaml(system, yamloptions)
            else:
                system = self.system_from_yaml(yamloptions)
        except KeyError:
            log.info("Config file doesn't contain any system info")
        
        if system:
            print(vars(system))
            return system
        else: 
            log.info("No system specified, quitting")
            sys.exit(1)

    def system_from_yaml(self, system_dict):
        """
        Super ugly code for making a system from the YAML file alone
        """
        #TODO Try to write better code. But first fix this thing.
        yamlSystem = YAMLSystem()
        for key, value in system_dict.iteritems(): 
            if key == "bins":
                print("TRYING TO SET %s FROM CONFIG FILE"%(key))
                setattr(init_system, "bin_mapper", bins_from_yaml_dict(value))
            #elif:
                # We want to parse target counts here
            #elif:
                # We want to parse any object here
            #elif:
                # We want to parse dtype here?
            else:
                # assume normal attribute
                try: 
                    print("TRYING TO SET %s FROM CONFIG FILE"%(key))
                    setattr(init_system, key, value)
                except:
                    log.info("Base system doesn't have propery %s"\
                          %(key))
                    pass
        # These are not parsed yet
        #self.pcoord_dtype = pcoord_dtype
        #self.bin_target_counts[...] = 10
        return yamlSystem


    def update_from_yaml(self, init_system, system_dict):
        """
        Updates the system built from the driver with the options
        from the YAML file.
        For now let's just overwrite anything.
        """
        #TODO Try to write better code. But first fix this thing.
        for key, value in system_dict.iteritems():
            if key == "bins":
                print("TRYING TO ADD BINS FROM CONFIG FILE")
                print("Bins of both were:")
                print(init_system.bin_mapper.boundaries)
                setattr(init_system, "bin_mapper", bins_from_yaml_dict(value))
                print(init_system.bin_mapper.boundaries)
            #elif:
                # We want to parse target counts here
            #elif:
                # We want to parse any object here
            #elif:
                # We want to parse dtype here?
            else:     
                try: 
                    print("TRYING TO SET %s FROM CONFIG FILE"%(key))
                    setattr(init_system, key, value)
                except:
                    log.info("Initial system doesn't have propery %s"\
                          %(key))
                    pass
        # These are not parsed yet
        #self.pcoord_dtype = pcoord_dtype
        #self.bin_target_counts[...] = 10
        return init_system

    def get_system_driver(self):
        if self._system is None:
            self._system = self.new_system_driver()
        return self._system
    
    def get_work_manager(self):
        return self.work_manager
    
    propagator = property(get_propagator)
    we_driver = property(get_we_driver)
    system = property(get_system_driver)
    data_manager = property(get_data_manager)
    sim_manager = property(get_sim_manager)
