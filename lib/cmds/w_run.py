from __future__ import division, print_function

import os, sys, logging, argparse, traceback
log = logging.getLogger('w_run')

import wemd

parser = wemd.rc.common_arg_parser(prog='w_run', description='start/continue a WEMD simulation')
parser.add_argument('--work-manager', dest='work_manager_name',
                    help='use the given work manager to distribute work among processors (serial, threads, processes, tcpip, zmq, '
                        +'or name a Python class; default: threads)')
mode_group = parser.add_mutually_exclusive_group()
mode_group.add_argument('--master', dest='mode', action='store_const', const='master',
                        help='Run as a WEMD master (responsible for coordinating WE and parallel propagation')
mode_group.add_argument('--node', dest='mode', action='store_const', const='node',
                        help='Run as a WEMD node coordinator (a communications broker between a remote master and local workers)')
mode_group.add_argument('--worker', dest='mode', action='store_const', const='worker',
                        help='Run as a WEMD worker (listening for work from a master)')
mode_group.set_defaults(mode='master') 
parser.add_argument('--oneseg', dest='only_one_segment', action='store_true',
                    help='only propagate one segment (useful for debugging propagators)')
parser.add_argument('--help-work-manager', dest='do_work_manager_help', action='store_true',
                    help='display help specific to the given work manager')

(args, aux_args) = parser.parse_known_args()

wemd.rc.config_logging(args, 'w_run')
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)

# Load the sim manager
sim_manager = wemd.rc.load_sim_manager(runtime_config)

# Load the work manager and have it parse remaining arguments
work_manager = sim_manager.load_work_manager()
aux_args = work_manager.parse_aux_args(aux_args, do_help=args.do_work_manager_help)
if aux_args:
    sys.stderr.write('unexpected command line argument(s) encountered: {}\n'.format(aux_args))
    sys.exit(os.EX_USAGE)

# Load remaining drivers
sim_manager.load_data_manager()
sim_manager.load_system_driver()
sim_manager.load_we_driver()
sim_manager.load_propagator()
sim_manager.load_plugins()

# Have the work manager perform any preparation prior to simulation (spawning clients, etc)
work_manager.startup()

# Have the sim manager perform an
log.debug('preparing run')
sim_manager.prepare_run()

try:
    log.debug('beginning run')
    sim_manager.run()
    
    log.debug('finalizing run')
    sim_manager.finalize_run()
finally:
    log.debug('back from sim manager')
    if not work_manager.shutdown_called:
        log.debug('work_manager.shutdown() not called -- calling directly')
        work_manager.shutdown()

if args.debug_mode:
    import threading, thread
    for thread in threading.enumerate():
        if 'MainThread' not in thread.name:
            log.debug('thread {!r} is still alive'.format(thread))
