from __future__ import print_function, division; __metaclass__ = type

import wemd

class WEMDTool:
    
    prog = None
    usage = None
    description = None
    epilog = None
    
    def __init__(self):
        self.wt2_config_required = False
    
    def add_args(self, parser):
        '''Add arguments specific to this tool to the given argparse parser.'''
        wemd.rc.add_args(parser)
    
    def process_args(self, args):
        '''Take argparse-processed arguments associated with this tool and deal
        with them appropriately (setting instance variables, etc)'''
        wemd.rc.process_args(args, config_required = self.wt2_config_required)
        
    def add_all_args(self, parser):
        '''Add arguments for all tools to the given parser.'''
        for cls in reversed(self.__class__.__mro__):
            try:
                fn = cls.__dict__['add_args']
            except KeyError:
                pass
            else:
                fn(self,parser)
    
    def process_all_args(self, args):
        '''Process arguments for all tools.'''
        for cls in reversed(self.__class__.__mro__):
            try:
                fn = cls.__dict__['process_args']
            except KeyError:
                pass
            else:
                fn(self,args)
                    
    def make_parser_and_process(self, prog=None, usage=None, description=None, epilog=None, args=None):
        '''A convenience function to create a parser, call add_all_args(), and then call process_all_args().
        The argument namespace is returned.'''
        import argparse
        prog = prog or self.prog
        usage = usage or self.usage
        description = description or self.description
        epilog = epilog or self.epilog
        parser = argparse.ArgumentParser(prog=prog, usage=usage, description=description, epilog=epilog)
        self.add_all_args(parser)
        args = parser.parse_args(args)
        self.process_all_args(args)
        return args
    
    def go(self):
        '''Perform the analysis associated with this object.'''
        raise NotImplementedError
    
    def main(self):
        '''A convenience function to make a parser, parse and process arguments, then call self.go()'''
        self.make_parser_and_process()
        self.go()
    
                
        