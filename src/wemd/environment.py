import os

ENV_STATE                = 'WE_STATE'
ENV_CURRENT_ITER         = 'WE_CURRENT_ITER'
ENV_CURRENT_SEG_ID       = 'WE_CURRENT_SEG_ID'
ENV_CURRENT_SEG_DATA_REF = 'WE_CURRENT_SEG_DATA_REF'
ENV_PARENT_ITER          = 'WE_PARENT_ITER'
ENV_PARENT_SEG_ID        = 'WE_PARENT_SEG_ID'
ENV_PARENT_SEG_DATA_REF  = 'WE_PARENT_SEG_DATA_REF'

EX_SUCCESS           = 0
EX_ERROR             = 1
EX_USAGE_ERROR       = 2
EX_ENVIRONMENT_ERROR = 3
EX_STATE_ERROR       = 4

class WEEnvironmentError(EnvironmentError):
    pass

class WEEnvironmentVarError(WEEnvironmentError):
    def __init__(self, e):
        self.message = 'WE environment variable not set: %s' % e

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

_data_manager = None
_current_segment = None

def get_data_manager():
    global _data_manager
    
    if _data_manager is not None: 
        return _data_manager
    
    source = get_we_variable(ENV_STATE)
    from we.data_managers import make_data_manager
    _data_manager = make_data_manager(source)
    _data_manager.restore_state()
    return _data_manager

def get_current_segment(load_parent = True):
    global _current_segment
    
    if _current_segment is not None:
        return _current_segment

    datamgr = get_data_manager()
    we_iter = long(get_we_variable(ENV_CURRENT_ITER))
    seg_id = long(get_we_variable(ENV_CURRENT_SEG_ID))
        
    _current_segment = get_data_manager().get_segment(we_iter, seg_id, 
                                                      load_parent)
    return _current_segment

__all__ = [name for name in dict(locals()) 
           if not name.startswith('_') 
           and name not in ('os',)]