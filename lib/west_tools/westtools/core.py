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


'''Core classes for creating WESTPA command-line tools'''


import sys, argparse, os
import westpa
import work_managers
from westpa import h5io

import logging
log = logging.getLogger(__name__)

class WESTToolComponent:
    '''Base class for WEST command line tools and components used in constructing tools'''
    
    def __init__(self):
        self.config_required = False
        self.include_args = {}
        self.arg_defaults = {}
        self.parser = None
        self.args = None
        
    def include_arg(self, argname):
        self.include_args[argname] = True
        
    def exclude_arg(self, argname):
        self.include_args[argname] = False
        
    def set_arg_default(self, argname, value):
        self.arg_defaults[argname] = value
        
    def add_args(self, parser):
        '''Add arguments specific to this component to the given argparse parser.'''
        pass
    
    def process_args(self, args):
        '''Take argparse-processed arguments associated with this component and deal
        with them appropriately (setting instance variables, etc)'''
        pass

    def add_all_args(self, parser):
        '''Add arguments for all components from which this class derives to the given parser,
        starting with the class highest up the inheritance chain (most distant ancestor).'''
        self.parser = parser
        for cls in reversed(self.__class__.__mro__):
            try:
                fn = cls.__dict__['add_args']
            except KeyError:
                pass
            else:
                fn(self,parser)
    
    def process_all_args(self, args):
        self.args = args
        '''Process arguments for all components from which this class derives,
        starting with the class highest up the inheritance chain (most distant ancestor).'''
        for cls in reversed(self.__class__.__mro__):
            try:
                fn = cls.__dict__['process_args']
            except KeyError:
                pass
            else:
                fn(self,args)

class WESTTool(WESTToolComponent):
    '''Base class for WEST command line tools'''
    
    prog = None
    usage = None
    description = None
    epilog = None
    
    def __init__(self):
        super(WESTTool,self).__init__()
                    
    def add_args(self, parser):
        '''Add arguments specific to this tool to the given argparse parser.'''
        westpa.rc.add_args(parser)
    
    def process_args(self, args):
        '''Take argparse-processed arguments associated with this tool and deal
        with them appropriately (setting instance variables, etc)'''
        westpa.rc.process_args(args, config_required = self.config_required)
                        
    def make_parser(self, prog=None, usage=None, description=None, epilog=None, args=None):
        prog = prog or self.prog
        usage = usage or self.usage
        description = description or self.description
        epilog = epilog or self.epilog
        parser = argparse.ArgumentParser(prog=prog, usage=usage, description=description, epilog=epilog,
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         conflict_handler='resolve')
        self.add_all_args(parser)
        return parser
            
    def make_parser_and_process(self, prog=None, usage=None, description=None, epilog=None, args=None):
        '''A convenience function to create a parser, call add_all_args(), and then call process_all_args().
        The argument namespace is returned.'''
        parser = self.make_parser(prog,usage,description,epilog,args)
        args = parser.parse_args(args)
        self.process_all_args(args)
        return args
    
    def go(self):
        '''Perform the analysis associated with this tool.'''
        raise NotImplementedError
    
    def main(self):
        '''A convenience function to make a parser, parse and process arguments, then call self.go()'''
        self.make_parser_and_process()
        self.go()
    
                
class WESTParallelTool(WESTTool):
    '''Base class for command-line tools parallelized with wwmgr. This automatically adds and processes
    wwmgr command-line arguments and creates a work manager at self.work_manager.'''

    def __init__(self, wm_env=None):
        super(WESTTool,self).__init__()
        self.work_manager = None
        self.wm_env = wm_env or work_managers.environment.default_env
        self.max_queue_len = None

    def make_parser_and_process(self, prog=None, usage=None, description=None, epilog=None, args=None):
        '''A convenience function to create a parser, call add_all_args(), and then call process_all_args().
        The argument namespace is returned.'''
        parser = self.make_parser(prog,usage,description,epilog,args)
        self.wm_env.add_wm_args(parser)
        
        args = parser.parse_args(args)
        self.wm_env.process_wm_args(args)
        
        # Instantiate work manager        
        self.work_manager = self.wm_env.make_work_manager()
        
        # Process args
        self.process_all_args(args)
        return args

    def add_args(self, parser):
        pgroup = parser.add_argument_group('parallelization options')
        pgroup.add_argument('--max-queue-length', type=int,
                            help='''Maximum number of tasks that can be queued. Useful to limit RAM use
                            for tasks that have very large requests/response. Default: no limit.''')
    
    def process_args(self, args):
        self.max_queue_len = args.max_queue_length
        log.debug('max queue length: {!r}'.format(self.max_queue_len))

    def go(self):
        '''Perform the analysis associated with this tool.'''
        raise NotImplementedError
    
    def main(self):
        '''A convenience function to make a parser, parse and process arguments, then run self.go() in the master process.'''
        self.make_parser_and_process()
        with self.work_manager:
            if self.work_manager.is_master:
                self.go()
            else:
                self.work_manager.run()

class WESTSubcommand(WESTToolComponent):
    '''Base class for command-line tool subcommands. A little sugar for making this 
    more uniform.'''
    
    subcommand = None
    help_text = None
    description = None
    
    def __init__(self, parent):
        self.parent = parent
        self.subparser = None
        
    def add_to_subparsers(self, subparsers):
        subparser = subparsers.add_parser(self.subcommand, help=self.help_text, description=self.description,
                                          formatter_class=argparse.RawDescriptionHelpFormatter,
                                          )
        self.add_all_args(subparser)
        subparser.set_defaults(west_subcommand=self)
        self.subparser = subparser

    def go(self):
        raise NotImplementedError

    @property
    def work_manager(self):
        '''The work manager for this tool. Raises AttributeError if this is not a parallel
        tool.'''
        return self.parent.work_manager

    
class _WESTSubcommandHelp(WESTSubcommand):
    subcommand = 'help'
    help_text = 'print help for this command or individual subcommands'
            
    def add_args(self, parser):
        parser.add_argument('command', nargs='?', choices=[subcommand.subcommand for subcommand in self.parent.subcommands])
    
    def process_args(self, args):
        self.command = args.command
        
    def go(self):
        if self.command is None:
            # Get parent help
            self.parent.parser.print_help(sys.stdout)
        else:
            self.parent._subcommand_instances[self.command].subparser.print_help(sys.stdout)
        sys.exit(0)
        
        

class WESTMasterCommand(WESTTool):
    '''Base class for command-line tools that employ subcommands'''
    
    subparsers_title = None
    subcommands = None
    
    include_help_command = True
    
    def __init__(self):
        super(WESTMasterCommand,self).__init__()
        self._subcommand = None
        self._subcommand_instances = {subcommand_class.subcommand: subcommand_class(self) for subcommand_class in self.subcommands}
        
        # Sanity checks
        if __debug__:
            for scclass in self.subcommands:
                assert scclass.subcommand, 'subcommand {!r} does not define a "subcommand" class variable'.format(scclass)
    
    def add_args(self, parser):
        subparsers = parser.add_subparsers(title=self.subparsers_title)
        if self.include_help_command:
            _WESTSubcommandHelp(self).add_to_subparsers(subparsers)
        for instance in self._subcommand_instances.values():
            instance.add_to_subparsers(subparsers)

    def process_args(self, args):
        try:
            self._subcommand = args.west_subcommand
        except AttributeError:
            # No subcommand given; display help
            print('Error: a command is required. See below.', file=sys.stderr)
            self.parser.print_help(sys.stderr)
            sys.exit(2)
        else:
            self._subcommand.process_all_args(args)
    
    def go(self):
        self._subcommand.go()
    
