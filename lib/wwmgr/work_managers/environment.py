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

'''Routines for configuring the work manager environment'''


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
    
    group_title = 'parallelization options'
    group_description = None
          
    env_prefix = 'WM'
    arg_prefix = 'wm'
    
    default_work_manager = 'serial'
    default_parallel_work_manager = 'processes'
    valid_work_managers = list(_available_work_managers.keys())
    
    def __init__(self, use_arg_prefixes=False, valid_work_managers=None):
        self.environ = os.environ
        self.args = None
        
        # copy from the class variable to permit modification
        # this can be modified to disable certain work managers, for whatever reason
        # mostly it's about having a valid list for enumeration
        self.valid_work_managers = valid_work_managers or list(self.valid_work_managers)
        
        self.use_arg_prefixes=use_arg_prefixes
        
    def env_name(self, name):
        return '{}_{}'.format(self.env_prefix, name.upper())
    
    def arg_name(self, name):
        if self.use_arg_prefixes:
            return '{}_{}'.format(self.arg_prefix, name)
        else:
            return name
    
    def arg_flag(self, name):
        if self.use_arg_prefixes:
            return '--{}-{}'.format(self.arg_prefix, re.sub('_','-',name))
        else:
            return '--{}'.format(re.sub('_','-',name))
            
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
        
        wm_group = parser.add_argument_group(self.group_title, self.group_description)
        wm_mutex = wm_group.add_mutually_exclusive_group()
        wm_mutex.add_argument(self.arg_flag('serial'), dest=self.arg_name('work_manager'), action='store_const', const='serial',
                              help='run in serial mode')
        wm_mutex.add_argument(self.arg_flag('parallel'), dest=self.arg_name('work_manager'), action='store_const',
                              const=self.default_parallel_work_manager, help='run in parallel mode (using {})'
                                                                             .format(self.default_parallel_work_manager))
        wm_mutex.add_argument(self.arg_flag('work_manager'), metavar='WORK_MANAGER',
                              choices=self.valid_work_managers,
                              help='''use the given work manager for parallel task distribution. Available
                                   work managers are {!r}; default is {!r}
                                   '''.format(tuple(self.valid_work_managers), self.default_work_manager))
            
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
    