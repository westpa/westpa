import os, sys
from optparse import OptionParser
import wemd, wemd.util.mpi
from wemd import Segment, WESimIter
from wemd.util.wetool import WECmdLineTool
from wemd.rc import EX_EXCEPTION_ERROR

import logging
log = logging.getLogger(__name__)

class WEMDRunTool(WECmdLineTool):
    def __init__(self):
        super(WEMDRunTool,self).__init__()
        
    def run(self):
        parser = OptionParser(usage = '%prog [options]')
        self.add_rc_option(parser)
        self.add_server_option(parser)
        (opts, args) = parser.parse_args()
        self.read_runtime_config(opts)
        
        sim_manager = self.load_sim_manager()
        
        if opts.hostname:
            try:
                sim_manager.set_hostname(opts.hostname)
            except Exception:
                sim_manager.shutdown(EX_EXCEPTION_ERROR)            
                log.error('unhandled exception in run() -- make sure worker is tcp server/client')
                sys.exit(EX_EXCEPTION_ERROR)    
                            
        sim_manager.load_sim_state()
        
        try:
            sim_manager.run()
        except Exception:
            log.error('unhandled exception in run()', exc_info = True)
            sim_manager.shutdown(EX_EXCEPTION_ERROR)
            sys.exit(EX_EXCEPTION_ERROR)
        else:
            log.info('WEMD run completed successfully')
            sim_manager.shutdown(0)
            sys.exit(0)
        finally:
            wemd.util.mpi.finalize_mpi()
        
if __name__ == '__main__':
    WEMDRunTool().run()

