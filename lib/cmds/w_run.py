from __future__ import division, print_function

import os, sys, logging, argparse, traceback, signal
log = logging.getLogger('w_run')

import wemd

parser = argparse.ArgumentParser('w_run', 'start/continue a WEMD simulation')
wemd.rc.add_args(parser)
parser.add_argument('--work-manager', dest='work_manager_name',
                    help='use the given work manager to distribute work among processors (serial, threads, processes, tcpip, zmq, '
                        +'or name a Python class; default: threads)')
parser.add_argument('--oneseg', dest='only_one_segment', action='store_true',
                    help='only propagate one segment (useful for debugging propagators)')
parser.add_argument('--kill-children', dest='do_kill_children', action='store_true',
                    help='''Forcibly kill child processes at program termination. (Risks data corruption; should only be used
                    for debugging.''')
parser.add_argument('--help-work-manager', dest='do_work_manager_help', action='store_true',
                    help='display help specific to the given work manager')

(args, aux_args) = parser.parse_known_args()
wemd.rc.process_args(args)

# Load the sim manager
sim_manager = wemd.rc.get_sim_manager()

# Load the work manager and have it parse remaining arguments
work_manager = wemd.rc.get_work_manager()
aux_args = work_manager.parse_aux_args(aux_args, do_help=args.do_work_manager_help)
if aux_args:
    sys.stderr.write('unexpected command line argument(s) encountered: {}\n'.format(aux_args))
    sys.exit(os.EX_USAGE)

# Load remaining drivers
system = wemd.rc.get_system_driver()
data_manager = wemd.rc.get_data_manager()
we_driver = wemd.rc.get_we_driver()
propagator = wemd.rc.get_propagator()

work_manager.propagator = propagator
propagator.system = system
data_manager.system = system
we_driver.system = system

sim_manager.work_manager = work_manager
sim_manager.data_manager = data_manager
sim_manager.system = system
sim_manager.propagator = propagator
sim_manager.we_driver = we_driver

sim_manager.load_plugins()

# Have the work manager perform any preparation prior to simulation (spawning clients, etc)
log.debug('calling work_manager.startup()')
work_manager.startup()

if work_manager.mode == 'master':
    log.debug('preparing run')
    sim_manager.prepare_run()
    
    try:
        log.debug('beginning run')
        sim_manager.run()
        
        log.debug('finalizing run')
        sim_manager.finalize_run()
    finally:
        log.debug('back from sim manager')
        if sys.exc_info != (None, None, None):
            log.error('unhandled exception; traceback follows')
            if not work_manager.shutdown_called:
                work_manager.shutdown(4)
        else:
            if not work_manager.shutdown_called:
                work_manager.shutdown(0)
                            
else:
    log.debug('running worker loop')
    work_manager.run_worker()

if wemd.rc.verbosity == 'debug':
    import threading, thread
    for thread in threading.enumerate():
        if 'MainThread' not in thread.name:
            log.debug('thread {!r} is still alive'.format(thread))

if args.do_kill_children:
    log.warning('killing children')
    prev_handler = signal.signal(signal.SIGTERM, signal.SIG_IGN)
    try:
        os.kill(-os.getpid(), signal.SIGTERM)
    finally:
        signal.signal(signal.SIGTERM, prev_handler)
