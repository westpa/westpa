# WEMD run control routines

import os
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