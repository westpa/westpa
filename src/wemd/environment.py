import os
from wemd.core.errors import WEEnvironmentError, WEEnvironmentVarError

ENV_RUNTIME_CONFIG       = 'WEMDRC'

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

_ex_names = dict((code, name) for (name, code) in locals().iteritems()
                 if name.startswith('EX_'))

def get_exit_code_name(code):
    return _ex_names.get(code, 'error %d' % code)

LOGGING_DEFAULT_FORMAT = 'WEMD %(levelname)-8s -- %(message)s'
LOGGING_DEFAULT_SETTINGS = {'': {'level': 'WARNING',
                                 'format': LOGGING_DEFAULT_FORMAT},
                            'wemd': {'level': 'WARNING'},
                            'sqlalchemy.engine': {'level': 'ERROR'}
                            }

def get_we_variable(varname, *args):
    if len(args) > 0:
        raise TypeError('unexpected positional argument encountered: %r' 
                        % args[1])
        
    try:
        return os.environ[varname]
    except KeyError:
        if args:
            return args[0]
        else:
            raise WEEnvironmentVarError(varname)

__all__ = [name for name in dict(locals()) 
           if not name.startswith('_') 
           and name not in ('os',)]