'''Routines for configuring the work manager environment'''

__metaclass__ = type

import os, re
from . import _available_work_managers

class WMEnvironment:
    '''A class to encapsulate the environment in which work managers are instantiated;
    this controls how environment variables and command-line arguments are used to
    set up work managers. This could be used to cleanly instantiate two work managers
    within one application, but is really more about providing facilities to make
    it easier for individual work managers to configure themselves according to
    precendence of configuration information:
      1. command-line arguments
      2. environment variables
      3. defaults
    '''
          
    env_prefix = 'WM'
    arg_prefix = 'wm'
    
    default_work_manager = 'processes'
    valid_work_managers = list(_available_work_managers.iterkeys())
    
    def __init__(self):
        self.environ = os.environ
        self.args = None
        
        # copy from the class variable
        # this can be modified to disable certain work managers, for whatever reason
        # mostly it's about having a valid list for enumeration
        self.valid_work_managers = list(self.valid_work_managers)
        
    def env_name(self, name):
        return '{}_{}'.format(self.env_prefix, name.upper())
    
    def arg_name(self, name):
        return '{}_{}'.format(self.arg_prefix, name)
    
    def arg_flag(self, name):
        return '--{}-{}'.format(self.arg_prefix, re.sub('_','-',name))
            
    def get_val(self, name, default=None, type_=None):
        envname = self.env_name(name)
        argname = self.arg_name(name)

        val = getattr(self.args,argname, None)
        if val is None:
            try:
                val = self.environ[envname]
            except KeyError:
                val = default
        
        if type_ is None:
            return val
        else:
            try:
                return type_(val)
            except ValueError as e:
                raise ValueError('cannot convert {!r} to {!r}: {!s}'.format(val,type_,e))
            
    def add_wm_args(self, parser):
        
        wm_group = parser.add_argument_group('work manager options')
        wm_group.add_argument(self.arg_flag('work_manager'), metavar='WORK_MANAGER',
                              choices=self.valid_work_managers,
                              help='use the given work manager for parallel task distribution (default: {})'
                                   .format(self.default_work_manager))
        wm_group.add_argument(self.arg_flag('n_workers'), metavar='N_WORKERS', type=int,
                              help='''Use up to N_WORKERS on this host, for work managers which support this option.
                                      Use 0 for a dedicated server. (Ignored by work managers which do not support
                                      this option.)''')
        
        for wm in self.valid_work_managers:
            _available_work_managers[wm].add_wm_args(parser,self)
            
    def process_wm_args(self, args):
        self.args = args
            
    def make_work_manager(self):
        '''Using cues from the environment, instantiate a pre-configured work manager.'''
        
        work_manager_name = self.get_val('work_manager','').lower()
        work_manager_name = work_manager_name or self.default_work_manager
        
        if work_manager_name not in self.valid_work_managers:
            raise ValueError('work manager {!r} is invalid or unavailable'.format(work_manager_name))
        else:
            return _available_work_managers[work_manager_name].from_environ(self)
        
default_env = WMEnvironment()
make_work_manager = default_env.make_work_manager
add_wm_args = default_env.add_wm_args
process_wm_args = default_env.process_wm_args
    