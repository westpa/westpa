from __future__ import division, print_function

import os, sys, logging, argparse, traceback
log = logging.getLogger('w_run')

import wemd
import work_managers
from work_managers import make_work_manager

parser = argparse.ArgumentParser('w_run', 'start/continue a WEMD simulation')
wemd.rc.add_args(parser)
parser.add_argument('--oneseg', dest='only_one_segment', action='store_true',
                    help='only propagate one segment (useful for debugging propagators)')

work_managers.environment.add_wm_args(parser)

args = parser.parse_args()
wemd.rc.process_args(args)
work_managers.environment.process_wm_args(args)
work_manager = make_work_manager()

# Load the sim manager and other drivers
sim_manager = wemd.rc.get_sim_manager(work_manager)
system = wemd.rc.get_system_driver()
data_manager = wemd.rc.get_data_manager()
we_driver = wemd.rc.get_we_driver()
propagator = wemd.rc.get_propagator()

propagator.system = system
data_manager.system = system
we_driver.system = system

sim_manager.data_manager = data_manager
sim_manager.system = system
sim_manager.propagator = propagator
sim_manager.we_driver = we_driver

with work_manager:
    if work_manager.is_master:
        work_manager.install_sigint_handler()
        sim_manager.load_plugins()
        
        log.debug('preparing run')
        sim_manager.prepare_run()
        
        try:
            log.debug('beginning run')
            sim_manager.run()
            
            log.debug('finalizing run')
            sim_manager.finalize_run()
        except KeyboardInterrupt:
            wemd.rc.pstatus('interrupted; shutting down')
        except:
            wemd.rc.pstatus('exception caught; shutting down')
            log.error(traceback.format_exc())
    else:
        work_manager.run()
