from __future__ import division, print_function

import os, sys, logging, argparse, traceback, signal
log = logging.getLogger('w_run')

import wemd

parser = argparse.ArgumentParser('w_run', 'start/continue a WEMD simulation')
wemd.rc.add_args(parser)
wemd.rc.add_work_manager_args(parser)
parser.add_argument('--oneseg', dest='only_one_segment', action='store_true',
                    help='only propagate one segment (useful for debugging propagators)')

(args, aux_args) = parser.parse_known_args()
wemd.rc.process_args(args, aux_args)

# Load the sim manager and other drivers
sim_manager = wemd.rc.get_sim_manager()
system = wemd.rc.get_system_driver()
data_manager = wemd.rc.get_data_manager()
we_driver = wemd.rc.get_we_driver()
propagator = wemd.rc.get_propagator()

propagator.system = system
data_manager.system = system
we_driver.system = system

sim_manager.work_manager = wemd.rc.get_work_manager()
sim_manager.data_manager = data_manager
sim_manager.system = system
sim_manager.propagator = propagator
sim_manager.we_driver = we_driver

work_manager = wemd.rc.get_work_manager()
mode = work_manager.startup()

sim_manager.load_plugins()

if work_manager.mode == work_manager.MODE_MASTER:
    
    log.debug('preparing run')
    sim_manager.prepare_run()
    
    try:
        log.debug('beginning run')
        sim_manager.run()
        
        log.debug('finalizing run')
        sim_manager.finalize_run()
        
        work_manager.shutdown(0)
    except:
        import traceback
        wemd.rc.pstatus('exception caught; shutting down')
        log.error(traceback.format_exc())
        work_manager.shutdown(4)

if wemd.rc.debug_mode:
    import threading, thread
    for thread in threading.enumerate():
        if 'MainThread' not in thread.name:
            log.debug('thread {!r} is still alive'.format(thread))
