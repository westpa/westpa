# WEMD run control routines

import os, sys, argparse
from wemd.util.config_dict import ConfigDict
from wemd.util import extloader

# Runtime config file management
ENV_RUNTIME_CONFIG  = 'WEMDRC'
RC_DEFAULT_FILENAME = 'wemd.cfg'

def read_config(filename = None):
    if filename is None:
        filename = RC_DEFAULT_FILENAME
    
    cdict = ConfigDict()
    cdict.read_config_file(filename)
    
    return cdict

def load_sim_manager(runtime_config):
    drivername = runtime_config.get('drivers.sim_manager', 'default')
    if drivername.lower() == 'default':
        from wemd.sim_manager import WESimManager
        return WESimManager(runtime_config)
    else:
        pathinfo = runtime_config.get_pathlist('drivers.module_path')
        return extloader.get_object(drivername,pathinfo)(runtime_config)
    
def common_arg_parser(prog=None,description=None):
    '''Return an argparse.ArgumentParser pre-loaded with --rcfile, --verbose, --debug,
    --profile, and --version'''
    import wemd
    parser = argparse.ArgumentParser(prog=prog,description=description)
    parser.add_argument('-r', '--rcfile', metavar='RCFILE', dest='run_config_file',
                        help='use RCFILE as the WEMD run-time configuration file (default: %s)' 
                              % RC_DEFAULT_FILENAME)
    parser.add_argument('--verbose', dest='verbose_mode', action='store_true',
                        help='emit extra information')
    parser.add_argument('--debug', dest='debug_mode', action='store_true',
                        help='enable extra checks and emit copious information')
    parser.add_argument('--profile', dest='profile_mode', action='store_true',
                        help='run this process under the Python profiler')
    parser.add_argument('--version', action='version', version='WEMD version %s' % wemd.version)
    return parser

def config_logging(args, tool_logger_name = None):
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
    
    if tool_logger_name:
        logging_config['loggers'][tool_logger_name] = {'handlers': ['console'], 'propagate': False}
        
    if args.debug_mode:
        logging_config['root']['level'] = 'DEBUG'
        logging_config['handlers']['console']['formatter'] = 'debug'
    elif args.verbose_mode:
        logging_config['root']['level'] = 'INFO'
    else:
        logging_config['root']['level'] = 'WARNING'
        
    
    logging.config.dictConfig(logging_config)
    logging_config['incremental'] = True

def run_optionally_profiled(func, args, kwargs, cmdline_args):
    '''Run func(*args, **kwargs).  If cmdline_args.profile_mode is True, run the same
    expression under the cProfile profiler.'''
    args = args or tuple()
    kwargs = kwargs or dict()
        
    if cmdline_args.profile_mode:
        import cProfile, pstats
        from wemd import log
        
        profile_file = '_wemd_profile_{:d}.dat'.format(os.getpid())
        log.info('writing profiling information to {}'.format(profile_file))
        
        gvars = globals()
        lvars = locals()
        lvars['func'] = func
        lvars['args'] = args
        lvars['kwargs'] = kwargs
        try:
            cProfile.runctx('func(*args, **kwargs)', gvars, lvars, profile_file)
        finally:
            stats = pstats.Stats(profile_file)
            stats.sort_stats('time')
            stats.print_stats()
    else:
        func(*args, **kwargs)
        
def default_cmdline_dispatch(func, args, kwargs, cmdline_args, log):
    '''Run func(*args, **kwargs) wrapped in the default WEMD exception handler, optionally with
    profiling (if cmdline_args.profile_mode is True).  log should be the top-level logger of
    the module responsible for calling func().  If an exception is caught, the program will
    terminate with exit code 1.'''

    args = args or tuple()
    kwargs = kwargs or dict()
    
    try:
        run_optionally_profiled(func, args, kwargs, cmdline_args)
    except KeyboardInterrupt:
        sys.stderr.write('Interrupted.\n')
        sys.exit(os.EX_TEMPFAIL)
    except Exception as e:
        import logging
        log.error(str(e))
        sys.stderr.write('ERROR: {!s}\n'.format(e))
        if log.isEnabledFor(logging.INFO):
            import traceback
            traceback.print_exc()
        sys.exit(1)
        
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

def get_exit_code_name(code):
    return _ex_names.get(code, 'error %d' % code)

__all__ = [name for name in dict(locals()) 
           if not name.startswith('_') 
           and name not in ('os',)]