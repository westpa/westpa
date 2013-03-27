from __future__ import print_function, division; __metaclass__ = type

import westpa
import work_managers

class WESTToolComponent:
    '''Base class for WEST command line tools and components used in constructing tools'''
    def __init__(self):
        self.config_required = False
        self.include_args = {}
        
    def include_arg(self, argname):
        self.include_args[argname] = True
        
    def exclude_arg(self, argname):
        self.include_args[argname] = False
        
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
        for cls in reversed(self.__class__.__mro__):
            try:
                fn = cls.__dict__['add_args']
            except KeyError:
                pass
            else:
                fn(self,parser)
    
    def process_all_args(self, args):
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
        import argparse
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

