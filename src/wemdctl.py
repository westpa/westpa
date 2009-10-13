import wemd
from wemd import Segment, WESimIter

from wemd.util.wetool import WECmdLineMultiTool
from wemd.environment import *

class WEMDCtlTool(WECmdLineMultiTool):
    def __init__(self):
        super(WEMDCtlTool,self).__init__()
        cop = self.command_parser
        
        cop.add_command('init', 'initialize a new WE simulation',
                        self.cmd_init, True)
        cop.add_command('run', 'run a WE simulation',
                        self.cmd_run, True)
        self.add_rc_option(cop)
        
    def cmd_init(self, args):
        parser = self.make_parser('SIM_CONFIG', 
                                  'Configure a new WEMD simulation according '
                                  +'to the directives in SIM_CONFIG')
        (opts, args) = parser.parse_args(args)
        if len(args) != 1:
            self.error_stream.write('a simulation configuration file is required\n')
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
        else:
            sim_config_file = args[0]
        
        self.connect_db()
        from wemd.core.we_sim import WESimDriver
        from wemd.util.config_dict import ConfigDict
        sim_config = ConfigDict()
        try:
            sim_config.read_config_file(sim_config_file)
        except IOError, e:
            self.log.debug('cannot open simulation config file', exc_info=True)
            self.error_stream.write('cannot open simulation config file: %s\n' % e)
            self.exit(EX_ENVIRONMENT_ERROR)
        
        sim_driver = WESimDriver(self.runtime_config)
        sim_driver.init_runtime()
        sim_driver.init_sim(sim_config)
        sim_driver.save_state()
        self.exit()
        
    def cmd_run(self, args):
        parser = self.make_parser(description = 'run the WE simulation')
        (opts, args) = parser.parse_args(args)
        
        
        
if __name__ == '__main__':
    WEMDCtlTool().run()

