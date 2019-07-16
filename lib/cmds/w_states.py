# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.



import sys, logging, numpy, argparse
import io
log = logging.getLogger('w_init')

import work_managers
from work_managers import make_work_manager

import westpa
from west import Segment
from west.states import BasisState, TargetState

EPS = numpy.finfo(numpy.float64).eps

parser = argparse.ArgumentParser('w_states', description='''\
Display or manipulate basis (initial) or target (recycling) states for a WEST simulation.  By default, states are
displayed (or dumped to files).  If ``--replace`` is specified, all basis/target states are replaced for the
next iteration.  If ``--append`` is specified, the given target state(s) are appended to the list for the
next iteration.

Appending basis states is not permitted, as this would require renormalizing basis state
probabilities in ways that may be error-prone. Instead, use ``w_states --show --bstate-file=bstates.txt``
and then edit the resulting ``bstates.txt`` file to include the new desired basis states, then use
``w_states --replace --bstate-file=bstates.txt`` to update the WEST HDF5 file appropriately. 
''')
westpa.rc.add_args(parser)
smgroup = parser.add_argument_group('modes of operation')
mode_group = smgroup.add_mutually_exclusive_group()
mode_group.add_argument('--show', dest='mode', action='store_const', const='show',
                        help='Display current basis/target states (or dump to files).')
mode_group.add_argument('--append', dest='mode', action='store_const', const='append',
                        help='Append the given basis/target states to those currently in use.')
mode_group.add_argument('--replace', dest='mode', action='store_const', const='replace',
                        help='Replace current basis/target states with those specified.')
parser.add_argument('--bstate-file', metavar='BSTATE_FILE',
                    help='''Read (--append/--replace) or write (--show) basis state names, probabilities, 
                    and data references from/to BSTATE_FILE.''')
parser.add_argument('--bstate', action='append', dest='bstates',
                    help='''Add the given basis state (specified as a string 'label,probability[,auxref]')
                    to the list of basis states (after those specified in --bstate-file, if any). This argument
                    may be specified more than once, in which case the given states are appended in the order
                    they are given on the command line.''')
parser.add_argument('--tstate-file', metavar='TSTATE_FILE',
                    help='''Read (--append/--replace) or write (--show) target state names 
                    and representative progress coordinates from/to TSTATE_FILE''')
parser.add_argument('--tstate', action='append', dest='tstates',
                    help='''Add the given target state (specified as a string 'label,pcoord0[,pcoord1[,...]]') to the
                    list of target states (after those specified in the file given by --tstates-from, if any).
                    This argument may be specified more than once, in which case the given states are appended
                    in the order they appear on the command line.''')
parser.set_defaults(mode='show')

work_managers.environment.add_wm_args(parser)
args = parser.parse_args()
westpa.rc.process_args(args)
work_managers.environment.process_wm_args(args)
work_manager = make_work_manager()

system = westpa.rc.get_system_driver()

with work_manager:
    if work_manager.is_master:
        data_manager = westpa.rc.get_data_manager()
        data_manager.open_backing(mode='a')
        sim_manager = westpa.rc.get_sim_manager()
        n_iter = data_manager.current_iteration
            
        assert args.mode in ('show', 'replace', 'append')
        if args.mode == 'show':
    
            basis_states = data_manager.get_basis_states(n_iter)
            if basis_states:
                bstate_file = sys.stdout if not args.bstate_file else open(args.bstate_file, 'wt')
                bstate_file.write('# Basis states for iteration {:d}\n'.format(n_iter))
                BasisState.states_to_file(basis_states, bstate_file)
            
            target_states = data_manager.get_target_states(n_iter)
            if target_states:
                tstate_file = sys.stdout if not args.tstate_file else open(args.tstate_file, 'wt')
                tstate_file.write('# Target states for iteration {:d}\n'.format(n_iter))
                TargetState.states_to_file(target_states, tstate_file)
                
        elif args.mode == 'replace':
            seg_index = data_manager.get_seg_index(n_iter)
            if (seg_index['status'] == Segment.SEG_STATUS_COMPLETE).any():
                print('Iteration {:d} has completed segments; applying new states to iteration {:d}'.format(n_iter,n_iter+1))
                n_iter += 1
            
            basis_states = []
            if args.bstate_file:
                basis_states.extend(BasisState.states_from_file(args.bstate_file))
            if args.bstates:
                for bstate_str in args.bstates:
                    fields = bstate_str.split(',')
                    label=fields[0]
                    probability=float(fields[1])
                    try:
                        auxref = fields[2]
                    except IndexError:
                        auxref = None
                    basis_states.append(BasisState(label=label,probability=probability,auxref=auxref))
                    
            if basis_states:
                # Check that the total probability of basis states adds to one
                tprob = sum(bstate.probability for bstate in basis_states)
                if abs(1.0 - tprob) > len(basis_states) * EPS:
                    pscale = 1/tprob
                    log.warning('Basis state probabilities do not add to unity; rescaling by {:g}'.format(pscale))
                    for bstate in basis_states:
                        bstate.probability *= pscale        
                
                # Assign progress coordinates to basis states
                sim_manager.get_bstate_pcoords(basis_states, n_iter)
                data_manager.create_ibstate_group(basis_states, n_iter)
                sim_manager.report_basis_states(basis_states)
                
            # Now handle target states
            target_states = []
            if args.tstate_file:
                target_states.extend(TargetState.states_from_file(args.tstate_file, system.pcoord_dtype))
            if args.tstates:
                tstates_strio = io.StringIO('\n'.join(args.tstates).replace(',', ' '))
                target_states.extend(TargetState.states_from_file(tstates_strio, system.pcoord_dtype))
                del tstates_strio
                
            if not target_states:
                westpa.rc.pstatus('No target states specified.')
            else:
                data_manager.save_target_states(target_states, n_iter)
                sim_manager.report_target_states(target_states)
                
            data_manager.update_iter_group_links(n_iter)
            
        else: # args.mode == 'append'
            if args.bstate_file or args.bstates:
                sys.stderr.write('refusing to append basis states; use --show followed by --replace instead\n')
                sys.exit(2)
            
            target_states = data_manager.get_target_states(n_iter)
            
            seg_index = data_manager.get_seg_index(n_iter)
            if (seg_index['status'] == Segment.SEG_STATUS_COMPLETE).any():
                print('Iteration {:d} has completed segments; applying new states to iteration {:d}'.format(n_iter,n_iter+1))
                n_iter += 1
            
            if args.tstate_file:
                target_states.extend(TargetState.states_from_file(args.tstate_file, system.pcoord_dtype))
            if args.tstates:
                tstates_strio = io.StringIO('\n'.join(args.tstates).replace(',', ' '))
                target_states.extend(TargetState.states_from_file(tstates_strio, system.pcoord_dtype))
                del tstates_strio
                
            if not target_states:
                westpa.rc.pstatus('No target states specified.')
            else:
                data_manager.save_target_states(target_states, n_iter)
                sim_manager.report_target_states(target_states)
                
            data_manager.update_iter_group_links(n_iter)
    else:
        work_manager.run()
