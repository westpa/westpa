from __future__ import division, print_function

import os, sys, logging, argparse, traceback
log = logging.getLogger('w_run')

import wemd
from work_managers import make_work_manager

#def remote_initialize(args):
#    import wemd
#    wemd.rc.process_args(args)
#    propagator = wemd.rc.get_propagator()
#    system = wemd.rc.get_system_driver()
#    propagator.system = system

parser = argparse.ArgumentParser('w_run', 'start/continue a WEMD simulation')
wemd.rc.add_args(parser)
parser.add_argument('--oneseg', dest='only_one_segment', action='store_true',
                    help='only propagate one segment (useful for debugging propagators)')

(args, aux_args) = parser.parse_known_args()
wemd.rc.process_args(args, aux_args)

# Load the sim manager and other drivers
work_manager = make_work_manager()
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

work_manager.startup()
try:
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
finally:
    work_manager.shutdown()
