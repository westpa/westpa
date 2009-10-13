from optparse import OptionParser
import wemd
from wemd import Segment, WESimIter

from wemd.backend_drivers import make_backend_driver
from wemd.work_managers import make_work_manager
from wemd.util.wetool import WECmdLineTool
from wemd.environment import *

class WEMDRunTool(WECmdLineTool):
    def __init__(self):
        super(WEMDRunTool,self).__init__()
        
    def run(self):
        parser = OptionParser(usage = '%prog [options]')
        self.add_rc_option(parser)
        (opts, args) = parser.parse_args()
        self.read_runtime_config(opts)
        for key in ('work_manager.driver', 'backend.driver'):
            self.runtime_config.require(key)
            
        backend_driver = make_backend_driver(self.runtime_config)
        backend_driver.initialize(self.runtime_config)
        work_manager = make_work_manager(self.runtime_config)
        work_manager.initialize(backend_driver, self.runtime_config)
        
        sim_driver = wemd.core.we_sim.WESimDriver(self.runtime_config)
        sim_driver.init_runtime()
        sim_driver.connect_db()
        sim_driver.restore_state()
        sim_driver.work_manager = work_manager
        work_manager.sim_driver = sim_driver 
        sim_driver.run_sim()
        self.exit(0)
        
        
if __name__ == '__main__':
    WEMDRunTool().run()

