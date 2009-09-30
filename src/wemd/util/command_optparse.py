import sys
from optparse import OptionParser, OptParseError

class NoCommandError(OptParseError):
    def __init__(self):
        OptParseError.__init__(self, 'no command given')

class InvalidCommandError(OptParseError):
    def __init__(self, command):
        OptParseError.__init__("no such command: '%s'" % command)

class Command(object):
    def __init__(self,
                 name, description, fn,
                 accepts_help_arg = False
                 ):
        
        self.name = name
        self.description = description
        self.accepts_help_arg = accepts_help_arg
        self.fn = fn
        
    def __call__(self, *args, **kwargs):
        return self.fn(*args, **kwargs)

class CommandOptionParser(OptionParser):
    invalid_command_exit_code = 2
    def __init__(self, *args, **kwargs):
        try:
            commands = kwargs.pop('commands')
        except KeyError:
            commands = []
        self.set_commands(commands)
        OptionParser.__init__(self, *args, **kwargs)
        self.disable_interspersed_args()
        
    def add_command(self, *args, **kwargs):
        try:
            cmdname = args[0].name
        except AttributeError:
            command = Command(*args, **kwargs)
        else:
            command = args[0]
            
        self.commands.append(command)
        self._command_index[command.name] = command
        
    def get_command(self, cmdname):
        return self._command_index[cmdname]
    
    def set_commands(self, commands):
        self.commands = []
        self._command_index = {}
        for command in commands:
            self.add_command(command)
        
    def format_command_help(self):
        import textwrap
        
        maxlen = max(len(cmd.name) for cmd in self.commands)
        wrapper = textwrap.TextWrapper(subsequent_indent = ' ' * (maxlen+3))
        return '\n'.join('%-*s - %s' 
                         % (maxlen, cmd.name, 
                            wrapper.fill(cmd.description))
                         for cmd in self.commands) + '\n'
            
    def format_help(self, formatter=None):
        txt = OptionParser.format_help(self, formatter)
        if self.commands:
            return '%s\ncommands:\n%s' % (txt, self.format_command_help())
        else:
            return txt
        
    def default_help_command(self, args = None, values = None):
        try:
            command = self._command_index[args[0]]
        except (KeyError,IndexError,TypeError):
            self.print_help()
        else:
            if command.accepts_help_arg:
                command(['--help'])
            else:
                self.print_help()
                
    def prescan(self, args = None):
        args = args or self._get_args()
        
        for arg in args:
            if arg.startswith('-'):
                continue
            else:
                break
        try:
            return self._command_index[arg]
        except KeyError:
            return None
        
    def parse_args(self, args = None, values = None):
        if 'help' not in self._command_index:
            self.add_command(Command('help', 'show help', 
                                     self.default_help_command,
                                     False))

        if values is None:
            values = self.get_default_values()
        else:
            values._update_loose(self.get_default_values().__dict__)
            
        (values, rargs) = OptionParser.parse_args(self, args, values)
        self.values = values 

        try:
            cmdname = rargs.pop(0)
        except IndexError:
            self.exit(self.invalid_command_exit_code,
                      'no command given\n%s' % self.format_help())
        
        try:
            self.command = self._command_index[cmdname]
        except KeyError:
            self.exit(self.invalid_command_exit_code,
                      "no such command: '%s'\n%s" 
                      % (cmdname, self.format_help()))            
        return (values, rargs)
