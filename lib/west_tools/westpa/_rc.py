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


import logging
log = logging.getLogger('westpa.rc')

import os, sys, errno, numpy, math, warnings

import westpa
from .yamlcfg import YAMLConfig
from .yamlcfg import YAMLSystem
from . import extloader
from work_managers import SerialWorkManager

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
        parsed_lists = boundary_lists[:]
        for iboundary, boundary in enumerate(boundary_lists):
            if boundary.__class__ == str:
                parsed_lists[iboundary] = parsePCV(boundary)[0]
            else: 
                parsed_lists[iboundary] = list(map((lambda x: float('inf') if (x if isinstance(x, str) else '').lower() == 'inf' else x), boundary))
        return mapper_type(parsed_lists)
    else:
        try:
            return mapper_type(**kwargs)
        except Exception:
            log.exception('exception instantiating mapper')
            raise

def parsePCV(pc_str):
    # Execute arbitrary code within a limited
    # scope to avoid nastyness. Stolen fully from 
    # other parts of the WESTPA code. 
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
    #return list(arr[...])
    return arr[...]

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
        self.config['args'] = {k:v for k,v in args.__dict__.items() if not k.startswith('_')}
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

        if self.verbosity == 'debug':
            warnings.resetwarnings()
            warnings.simplefilter("default")
            logging.captureWarnings(True)
        else:
            if not sys.warnoptions:
                warnings.simplefilter("ignore")
            logging.captureWarnings(False)
        
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
        Returns a new system object either from the driver OR from the YAML
        file. Currently builds off of system then updates with YAML
        overwriting the previous settings if both are specified.

        There are no default settings, all settings MUST be specified
        in the config YAML file OR the system driver. 

        Settings that are specified here are: 
          Progress coordinate settings:
            Progress coordinate dimensionality
            Progress coordinate length/number of data points in each tau
            Progress coordinate data type 
          Bin settings: 
            Bin mapper type
            Bins 
            Target simulation counts for each bin 
          Generic setting for flexibility:
            Both paths allow for generic attribute setting for
            future flexibility.
        '''

        # First step is to see if we have a driver specified
        # if so load that one first, setting the defaults for
        # YAML to possibly modify. 

        # I will keep this as is for now since I like the idea
        # of being able to set some defaults from the system 
        # driver and then just modify the YAML instead for basic
        # stuff. This might make it easier to set up systems since
        # playing around with the YAML format is easier. 
        
        # Still still has the end-user issue of possibly confusing 
        # people, maybe KISS is better, not sure. I'll keep this 
        # for development purposes if nothing else.

        system = None
        # Get method checks for us  
        sysdrivername = self.config.get(['west', 'system', 'driver']) 
        if not sysdrivername:
            # Warn user that driver is not specified
            log.info("System driver not specified")
        else:
            log.info('loading system driver %r' % sysdrivername)
            system = extloader.get_object(sysdrivername)(rc=self)
            log.debug('loaded system driver {!r}'.format(system))        
            system.initialize()
        # Second let's see if we have info in the YAML file 
        yamloptions = self.config.get(['west','system','system_options'])
        if not yamloptions:
            # Same here, yaml file doesn't have sys info
            log.info("Config file doesn't contain any system info")
        else:
             log.info("Loading system options from configuration file")
             if system:
                 system = self.update_from_yaml(system, yamloptions)
             else:
                 system = self.system_from_yaml(yamloptions)
        
        if system:
            if not yamloptions:
                print("System is being built only off of the system driver")
            return system
        else: 
            log.info("No system specified! Exiting program.")
            # Gracefully exit
            raise ValueError("No system defined!")

    def system_from_yaml(self, system_dict):
        """
        System builder directly from the config YAML file. 

        Arguments: 
          system_dict (dict): Parsed YAML file as a dictionary, parsed by 
            PyYAML by default.

        Returns: 
          A modified WESTSystem object as defined in yamlcfg.py with 
          the parsed settings from the config file.
        """
        
        yamlSystem = YAMLSystem()
        print("System building only off of the configuration file")  
        # Now for the building of the system from YAML we need to use
        # require for these settings since they are musts. 
        
        # First basic pcoord settings
        ndim  = self.config.require(\
            ['west', 'system', 'system_options', 'pcoord_ndim'])
        plen  = self.config.require(\
            ['west', 'system', 'system_options', 'pcoord_len'])
        # Dtype needs to be ran as code from YAML file, document YAML code execution syntax
        # somewhere
        ptype = self.config.require(\
            ['west', 'system', 'system_options', 'pcoord_dtype'])
        # Bins
        bins_obj = self.config.require(\
            ['west', 'system', 'system_options', 'bins'])
        trgt_cnt = self.config.require(\
            ['west', 'system', 'system_options', 'bin_target_counts'])
        # Now add the parsed settings to the system
        mapper   = bins_from_yaml_dict(bins_obj)
        setattr(yamlSystem, 'pcoord_ndim', ndim)
        setattr(yamlSystem, 'pcoord_len', plen)
        setattr(yamlSystem, 'pcoord_dtype', ptype)
        setattr(yamlSystem, 'bin_mapper', mapper)
        # Check if the supplied target count object is 
        # an iterable or not, 

        # This one I designed in a way that you can either
        # directly supply an iterable that has the correct
        # size OR an integer. Anything else will fail. 
        
        # Main issue: For higher bin dimensions this 
        # is tough since the bins might not correspond
        # properly to the flattened array you are supplying.

        # I might just scrap the iterable later and only allow
        # integers.
        if hasattr(trgt_cnt, "__iter__"):
            assert len(trgt_cnt) == mapper.nbins, \
              "Count iterable size doesn't match the number of bins"
            trgt_cnt_arr = trgt_cnt
        else:
            assert trgt_cnt == numpy.int(trgt_cnt), \
              "Counts are not integer valued, ambiguous input"
            trgt_cnt_arr    = numpy.zeros(mapper.nbins)
            trgt_cnt_arr[:] = trgt_cnt
        setattr(yamlSystem, 'bin_target_counts', trgt_cnt_arr)

        # Attach generic attribute to system 
        for attr in system_dict.keys():
            if not hasattr(yamlSystem, attr):
                setattr(yamlSystem, attr, system_dict[attr])

        # Return complete system
        return yamlSystem


    def update_from_yaml(self, init_system, system_dict):
        """
        Updates the system built from the driver with the options
        from the YAML file. For now it overwrites everything specified
        in the system driver.

        Arguments: 
          system_dict (dict): Parsed YAML file as a dictionary, parsed by 
            PyYAML by default.
          init_system (WESTSystem): System returned by the driver.

        Returns: 
          A modified WESTSystem object with settings from the system
          driver and the config YAML file.
        """
        
        # First we want to overwrite whatever we have from the YAML
        # file.
        print("Updating system with the options from the configuration file")
        for key, value in system_dict.items():
            if key == 'pcoord_ndim':
                self.overwrite_option(init_system, key, value)
            elif key == 'pcoord_len':
                self.overwrite_option(init_system, key, value)
            elif key == 'pcoord_dtype':
                self.overwrite_option(init_system, key, value)
            elif key == "bins":
                self.overwrite_option(init_system, 'bin_mapper',\
                     bins_from_yaml_dict(value))
        # Target counts have to be parsed after we have a mapper in
        # place
        try: 
            trgt_cnt = system_dict['bin_target_counts']
            if hasattr(trgt_cnt, "__iter__"):
                assert len(trgt_cnt) == init_system.bin_mapper.nbins, \
                  "Count iterable size doesn't match the number of bins"
                trgt_cnt_arr = trgt_cnt
            else:
                assert trgt_cnt == numpy.int(trgt_cnt), \
                  "Counts are not integer valued, ambiguous input"
                trgt_cnt_arr    = numpy.zeros(init_system.bin_mapper.nbins)
                trgt_cnt_arr[:] = int(trgt_cnt)
            self.overwrite_option(init_system, 'bin_target_counts', trgt_cnt_arr)
        except KeyError:
             pass
        # The generic attribute settings added here
        for attr in system_dict.keys():
            if not hasattr(init_system, attr):
                setattr(init_system, attr, system_dict[attr])
        return init_system

    def overwrite_option(self, system, key, value):
        if hasattr(system, key):
            log.info("Overwriting system option: %s"%key)
        setattr(system, key, value)

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
