import argparse
import logging
import traceback

import westpa
import westpa.work_managers as work_managers
from westpa.work_managers import make_work_manager


log = logging.getLogger('w_run')


def entry_point():
    parser = argparse.ArgumentParser('w_run', 'start/continue a WEST simulation')
    westpa.rc.add_args(parser)
    parser.add_argument(
        '--oneseg',
        dest='only_one_segment',
        action='store_true',
        help='only propagate one segment (useful for debugging propagators)',
    )

    work_managers.environment.add_wm_args(parser)

    args = parser.parse_args()
    westpa.rc.process_args(args)
    work_managers.environment.process_wm_args(args)
    work_manager = westpa.rc.work_manager = make_work_manager()

    # Load the sim manager and other drivers
    sim_manager = westpa.rc.get_sim_manager()
    system = westpa.rc.get_system_driver()
    data_manager = westpa.rc.get_data_manager()
    we_driver = westpa.rc.get_we_driver()
    propagator = westpa.rc.get_propagator()

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
                westpa.rc.pstatus('interrupted; shutting down')
            except Exception:
                westpa.rc.pstatus('exception caught; shutting down')
                log.error(traceback.format_exc())
        else:
            work_manager.run()


if __name__ == '__main__':
    entry_point()
