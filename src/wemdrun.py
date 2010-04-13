import os, sys
from optparse import OptionParser
import wemd, wemd.util.mpi
from wemd import Segment, WESimIter
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
        for key in ('backend.driver',):
            self.runtime_config.require(key)
            
        sim_manager = wemd.sim_managers.make_sim_manager(self.runtime_config)    
        sim_manager.restore_state()
        
        try:
            sim_manager.run()
        except:
            wemd.util.mpi.abort_mpi(1)
        finally:
            wemd.util.mpi.finalize_mpi()
        
if __name__ == '__main__':
    WEMDRunTool().run()

