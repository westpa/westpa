import argparse
import io
import logging
import sys

import numpy as np

import westpa

from westpa.core.states import BasisState, TargetState

import westpa.work_managers as work_managers
from westpa.work_managers import make_work_manager


log = logging.getLogger('w_init')
EPS = np.finfo(np.float64).eps


def entry_point():

    parser = argparse.ArgumentParser(
        'w_init',
        description='''\
    Initialize a new WEST simulation, creating the WEST HDF5 file and preparing the first
    iteration's segments.

    Initial states are generated from one or more "basis states" which are specified either in a file specified with --bstates-from, or
    by one or more "--bstate" arguments. If neither --bstates-from nor at least one "--bstate" argument is provided, then a default
    basis state of probability one identified by the state ID zero and label "basis" will be created (a warning will be printed in this
    case, to remind you of this behavior, in case it is not what you wanted).

    Target states for (non-equilibrium) steady-state simulations are specified either in a file specified with --tstates-from, or
    by one or more --tstate arguments. If neither --tstates-from nor at least one --tstate argument is provided, then an equilibrium
    simulation (without any sinks) will be performed.
    ''',
    )
    westpa.rc.add_args(parser)
    parser.add_argument('--force', dest='force', action='store_true', help='Overwrite any existing simulation data')
    parser.add_argument(
        '--bstate-file',
        '--bstates-from',
        metavar='BSTATE_FILE',
        help='Read basis state names, probabilities, and (optionally) data references from BSTATE_FILE.',
    )
    parser.add_argument(
        '--bstate',
        action='append',
        dest='bstates',
        help='''Add the given basis state (specified as a string 'label,probability[,auxref]')
                        to the list of basis states (after those specified in --bstates-from, if any). This argument
                        may be specified more than once, in which case the given states are appended in the order
                        they are given on the command line.''',
    )
    parser.add_argument(
        '--tstate-file',
        '--tstates-from',
        metavar='TSTATE_FILE',
        help='''Read target state names and representative progress coordinates from TSTATE_FILE. WESTPA uses the
                        representative progress coordinate of a target state and converts the **entire** bin
                        containing that progress coordinate into a recycling sink.''',
    )
    parser.add_argument(
        '--tstate',
        action='append',
        dest='tstates',
        help='''Add the given target state (specified as a string 'label,pcoord0[,pcoord1[,...]]') to the
                        list of target states (after those specified in the file given by --tstates-from, if any).
                        This argument may be specified more than once, in which case the given states are appended
                        in the order they appear on the command line.''',
    )

    parser.add_argument(
        '--sstate-file',
        '--sstates-from',
        metavar='SSTATE_FILE',
        help='Read start state names, probabilities, and (optionally) data references from SSTATE_FILE.',
    )
    parser.add_argument(
        '--sstate',
        action='append',
        dest='sstates',
        help='''Add the given start state (specified as a string 'label,probability[,auxref]')
                        to the list of start states (after those specified in --sstates-from, if any). This argument
                        may be specified more than once, in which case the given states are appended in the order
                        they are given on the command line.''',
    )

    parser.add_argument(
        '--segs-per-state',
        type=int,
        metavar='N',
        default=1,
        help='''Initialize N segments from each basis state (default: %(default)s).''',
    )
    parser.add_argument(
        '--no-we',
        '--shotgun',
        dest='shotgun',
        action='store_true',
        help='''Do not run the weighted ensemble bin/split/merge algorithm on newly-created segments.''',
    )

    # TODO: Does this belong here or not? I like that it's parsing arguments, which is the purpose of entry_point.
    #   I don't necessarily like that it's setting state across different parts of the program.

    work_managers.environment.add_wm_args(parser)
    args = parser.parse_args()
    westpa.rc.process_args(args)
    work_managers.environment.process_wm_args(args)

    initialize(
        args.tstates,
        args.tstate_file,
        args.bstates,
        args.bstate_file,
        args.sstates,
        args.sstate_file,
        args.segs_per_state,
        args.shotgun,
    )


def initialize(tstates, tstate_file, bstates, bstate_file, sstates=None, sstate_file=None, segs_per_state=1, shotgun=False):
    """
    Initialize a WESTPA simulation.

    tstates : list of str

    tstate_file : str

    bstates : list of str

    bstate_file : str

    sstates : list of str

    sstate_file : str

    segs_per_state : int

    shotgun : bool
    """

    westpa.rc.work_manager = work_manager = make_work_manager()

    system = westpa.rc.get_system_driver()
    sim_manager = westpa.rc.get_sim_manager()

    data_manager = westpa.rc.get_data_manager()

    data_manager.system = system

    with work_manager:
        if work_manager.is_master:
            # Process target states
            target_states = []
            if tstate_file:
                target_states.extend(TargetState.states_from_file(tstate_file, system.pcoord_dtype))
            if tstates:
                tstates_strio = io.StringIO('\n'.join(tstates).replace(',', ' '))
                target_states.extend(TargetState.states_from_file(tstates_strio, system.pcoord_dtype))
                del tstates_strio

            # Process basis states
            basis_states = []
            if bstate_file:
                basis_states.extend(BasisState.states_from_file(bstate_file))
            if bstates:
                for bstate_str in bstates:
                    fields = bstate_str.split(',')
                    label = fields[0]
                    probability = float(fields[1])
                    try:
                        auxref = fields[2]
                    except IndexError:
                        auxref = None
                    basis_states.append(BasisState(label=label, probability=probability, auxref=auxref))

            # Process the list of start states, creating a BasisState from each
            start_states = []
            if sstate_file:
                start_states.extend(BasisState.states_from_file(sstate_file))
            if sstates:
                for sstate_str in sstates:
                    fields = sstate_str.split(',')
                    label = fields[0]
                    probability = float(fields[1])
                    try:
                        auxref = fields[2]
                    except IndexError:
                        auxref = None
                    start_states.append(BasisState(label=label, probability=probability, auxref=auxref))

            if not basis_states:
                log.error('At least one basis state is required')
                sys.exit(3)

            # Check that the total probability of basis states adds to one
            bstate_prob, sstate_prob = (
                sum(bstate.probability for bstate in basis_states),
                sum(sstate.probability for sstate in start_states),
            )
            # tprob = sum(bstate.probability for bstate in basis_states)
            tprob = bstate_prob + sstate_prob
            if abs(1.0 - tprob) > len(basis_states) * EPS:
                pscale = 1 / tprob
                log.warning(
                    'Basis state probabilities do not add to unity (basis: {:.2f}, start states: {:.2f}); rescaling by {:g}. If using start states, some rescaling is normal.'.format(
                        bstate_prob, sstate_prob, pscale
                    )
                )
                for bstate in basis_states:
                    bstate.probability *= pscale
                for sstate in start_states:
                    sstate.probability *= pscale

            # Prepare simulation
            sim_manager.initialize_simulation(
                basis_states=basis_states,
                target_states=target_states,
                start_states=start_states,
                segs_per_state=segs_per_state,
                suppress_we=shotgun,
            )
        else:
            work_manager.run()


if __name__ == '__main__':
    entry_point()
