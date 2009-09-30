import sys, os, re, logging
from optparse import OptionParser
from command_optparse import CommandOptionParser
from wemd.environment import *

class WECmdLineTool(object):
    state_short_option = '-S'
    state_long_option  = '--state'
    seg_id_short_option = '-s'
    seg_id_long_option = '--seg-id'
    we_iter_short_option = '-i'
    we_iter_long_option  = '--iter'
    
    usage = '%prog [options] COMMAND [...]'
    description = None
    
    def __init__(self):
                
        self.input_stream  = sys.stdin
        self.output_stream = sys.stdout
        self.error_stream  = sys.stderr
        
        self.datamgr = None
        self.current_iteration = None
        self.current_segment = None
        
        logging.basicConfig(level=logging.WARN)
        
    def exit(self, code=0):
        sys.exit(code)
                
    def add_state_option(self, parser):
        parser.add_option(self.state_short_option, self.state_long_option, 
                          dest='state',
                          help='simulation status/data is located in STATE '
                              +'(default: taken from $%s)' % ENV_STATE)
        
    def add_seg_id_option(self, parser):
        parser.add_option(self.seg_id_short_option, self.seg_id_long_option,
                          dest='seg_id', type='long',
                          help='use SEG_ID as the applicable segment (default: '
                              +'taken from $%s' % ENV_CURRENT_SEG_ID)
        
    def add_we_iter_option(self, parser, source='state'):
        if source == 'environment':
            src_msg = '$' + ENV_CURRENT_ITER
        else:
            src_msg = 'stored state'
            
        parser.add_option(self.we_iter_short_option, self.we_iter_long_option,
                          dest='we_iter', type='long',
                          help='use WE_ITER as the applicable iteration '
                              +'(default: taken from %s)' % src_msg)
        
    def format_env_error_message(self, msg, optname, envname):
        return (msg.rstrip() + ' not specified on command line (%s) ' % optname
                + 'or in $%s' % envname) 
        
    def connect_datamgr(self, opts = None):
        from we.data_managers import make_data_manager
        try:
            source = opts.state
        except AttributeError:
            source = None
            
        if source is None:
            try:
                source = get_we_variable(ENV_STATE)
            except WEEnvironmentVarError:
                self.error_stream.write(
                    self.format_env_error_message('state/data source',
                                                  self.state_long_option,
                                                  ENV_STATE) + '\n')
                self.exit(EX_ENVIRONMENT_ERROR)
        self.datamgr = make_data_manager(source)
            
    def load_state(self, opts = None):
        if not self.datamgr:
            self.connect_datamgr(opts)
        self.datamgr.restore_state()
        
    def set_current_iteration(self, opts = None):
        try:
            we_iter = opts.we_iter
        except AttributeError:
            we_iter = None
            
        if we_iter is None:
            try:
                we_iter = long(get_we_variable(ENV_CURRENT_ITER))
            except WEEnvironmentVarError:
                try:
                    we_iter = self.datamgr.current_iteration
                except AttributeError:
                    self.error_stream.write(
                        self.format_env_error_message('current iteration', 
                                                      self.we_iter_long_option,
                                                      ENV_CURRENT_ITER) + '\n')
                    self.exit(EX_ENVIRONMENT_ERROR)
        self.current_iteration = we_iter
        
    def load_current_segment(self, opts = None):
        self.set_current_iteration(opts)
        try:
            seg_id = opts.seg_id
        except AttributeError:
            pass
        
        if seg_id is None:
            try:
                seg_id = long(get_we_variable(ENV_CURRENT_SEG_ID))
            except WEEnvironmentError:
                self.error_stream.write(
                    self.format_env_error_message('current segment ID', 
                                                  self.seg_id_long_option,
                                                  ENV_CURRENT_SEG_ID) + '\n')
                self.exit(EX_ENVIRONMENT_ERROR)
                
        self.current_segment = self.datamgr.get_segment(self.current_iteration,
                                                        seg_id)

    def run(self):
        raise NotImplementedError

class WECmdLineMultiTool(WECmdLineTool):
    
    def __init__(self):
        super(WECmdLineMultiTool,self).__init__()
        self.command_parser = CommandOptionParser(usage = self.usage,
                                                  description = self.description)
        
        self.global_opts = None
        self.command = None
        self.command_args = None

    def make_parser(self, usage_tail = '[options]', description = None):
        cmdname = self.command.name
        description = description or self.command.description or None
        parser = OptionParser(usage='%%prog %s %s' % (self.command.name,
                                                      usage_tail),
                              description = description)
        return parser

    def cmd_help(self, args):
        try:
            command = self.command_parser.get_command(args[0])
        except (KeyError,IndexError,TypeError):
            self.command_parser.print_help()
        else:
            if command.accepts_help_arg:
                self.command = command
                command(['--help'])
            else:
                self.command_parser.print_help()
        self.exit(0)
        
    def run(self):
        self.command_parser.add_command('help', 'show help', self.cmd_help, 
                                        False)
        (self.global_opts, self.command_args) = self.command_parser.parse_args()
        self.command = self.command_parser.command
        self.command(self.command_args)
