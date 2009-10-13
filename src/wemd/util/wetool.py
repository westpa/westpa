import sys, os, re, logging
from optparse import OptionParser
from command_optparse import CommandOptionParser
from wemd.environment import *
from wemd.util.config_dict import ConfigDict
from wemd.util.miscfn import logging_level_by_name


class WECmdLineTool(object):
    RC_SHORT_OPTION = '-r'
    RC_LONG_OPTION  = '--rcfile'
    RC_DEFAULT      = 'run.cfg'

    @classmethod
    def add_rc_option(cls, parser):
        parser.add_option(cls.RC_SHORT_OPTION, cls.RC_LONG_OPTION,
                          dest='rcfile', 
                          help='runtime configuration is in RCFILE '
                              +'(default: %s)' % cls.RC_DEFAULT,
                          default=cls.RC_DEFAULT)
        
    usage = '%prog [options] COMMAND [...]'
    description = None
    
    def __init__(self):    
        self.input_stream  = sys.stdin
        self.output_stream = sys.stdout
        self.error_stream  = sys.stderr
        self.runtime_config = None
        self.dbengine = None
        self.DBSession = None
        self.log = logging.getLogger('%s.%s' 
                                     % (__name__, self.__class__.__name__))
        
    def read_runtime_config(self, opts):
        self.runtime_config = ConfigDict()
        try:
            self.runtime_config.read_config_file(opts.rcfile)
        except IOError, e:
            self.error_stream.write('cannot open runtime config file: %s\n' % e)
            self.exit(EX_RTC_ERROR)
        self.configure_logging()
        
    def connect_db(self):
        self.runtime_config.require('data.db.url')
        db_url = self.runtime_config['data.db.url']
        self.log.info('connecting to %r' % db_url)
        
        import wemd.data_manager
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        
        self.dbengine = create_engine(db_url)
        self.DBSession = sessionmaker(bind=self.dbengine,
                                      autocommit = True,
                                      autoflush = False)
        self.runtime_config['data.db.engine'] = self.dbengine
        self.runtime_config['data.db.sessionmaker'] = self.DBSession
        
    def configure_logging(self):
        #logging.basicConfig()
        logkeys = set(k for k in self.runtime_config
                      if k.startswith('logging.'))
        
        logger_settings = LOGGING_DEFAULT_SETTINGS.copy()
        for key in logkeys:
            item = key[8:]
            value = self.runtime_config[key]
            try:
                (logger_name, param) = item.rsplit('.', 1)
            except ValueError:
                continue
            else:
                if param not in ('level', 'format', 'handler'):
                    continue
            if logger_name == 'root': logger_name = ''
            try:
                logger_settings[logger_name][param] = value
            except KeyError:
                logger_settings[logger_name] = {param: value}
                
        for (logger_name, settings) in logger_settings.iteritems():                
            level = logging_level_by_name(settings.get('level', 'NOTSET'))
            logger = logging.getLogger(logger_name)
            logger.setLevel(level)
            try:
                format = settings['format']
            except KeyError:
                pass
            else:
                handler = logging.StreamHandler()
                handler.setFormatter(logging.Formatter(format))
                logger.addHandler(handler)
                logger.propagate = False
                            
    def exit(self, code=0):
        self.log.info('exiting with code %d (%s)' 
                       % (code, get_exit_code_name(code)))
        sys.exit(code)

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
        self.read_runtime_config(self.global_opts)
        self.command = self.command_parser.command
        self.command(self.command_args)
