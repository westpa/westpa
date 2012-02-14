from __future__ import division, print_function

import os, sys, logging, numpy, operator, argparse
import cStringIO
from itertools import izip
log = logging.getLogger('w_init')

import wemd
from wemd.segment import Segment
from wemd.work_managers import ops
from wemd.states import BasisState, TargetState, InitialState, pare_basis_initial_states

EPS = numpy.finfo(numpy.float64).eps

parser = argparse.ArgumentParser('w_init', description='''\
Initialize a new WEMD simulation, creating the WEMD HDF5 file and preparing the first
iteration's segments.

Initial states are generated from one or more "basis states" which are specified either in a file specified with --bstates-from, or
by one or more "--bstate" arguments. If neither --bstates-from nor at least one "--bstate" argument is provided, then a default
basis state of probability one identified by the state ID zero and label "basis" will be created (a warning will be printed in this
case, to remind you of this behavior, in case it is not what you wanted).

Target states for (non-equilibrium) steady-state simulations are specified either in a file specified with --tstates-from, or
by one or more --tstate arguments. If neither --tstates-from nor at least one --tstate argument is provided, then an equilibrium
simulation (without any sinks) will be performed. 
''')
wemd.rc.add_args(parser)
wemd.rc.add_work_manager_args(parser)
parser.add_argument('--force', dest='force', action='store_true',
                         help='Overwrite any existing simulation data')
parser.add_argument('--bstates-from', metavar='BSTATE_FILE',
                    help='Read basis state names, probabilities, and (optionally) data references from BSTATE_FILE.')
parser.add_argument('--bstate', action='append', dest='bstates',
                    help='''Add the given basis state (specified as a string 'label,probability[,auxref]')
                    to the list of basis states (after those specified in --bstates-from, if any). This argument
                    may be specified more than once, in which case the given states are appended in the order
                    they are given on the command line.''')
parser.add_argument('--tstates-from', metavar='TSTATE_FILE',
                    help='Read target state names and representative progress coordinates from TSTATE_FILE')
parser.add_argument('--tstate', action='append', dest='tstates',
                    help='''Add the given target state (specified as a string 'label,pcoord0[,pcoord1[,...]]') to the
                    list of target states (after those specified in the file given by --tstates-from, if any).
                    This argument may be specified more than once, in which case the given states are appended
                    in the order they appear on the command line.''')
parser.add_argument('--segs-per-state', type=int, metavar='N', default=1,
                    help='''Initialize N segments from each basis state (default: %(default)s).''')
parser.add_argument('--run-we', action='store_true',
                    help='''Run the weighted ensemble bin/split/merge algorithm on newly-created segments.''')

(args, aux_args) = parser.parse_known_args()
wemd.rc.process_args(args, aux_args)
system = wemd.rc.get_system_driver()
sim_manager = wemd.rc.get_sim_manager()
propagator = wemd.rc.get_propagator()
h5file = wemd.rc.config.get_path('data.h5file')
data_manager = wemd.rc.get_data_manager()
data_manager.we_h5filename = h5file
data_manager.system = system
we_driver = wemd.rc.get_we_driver()

work_manager = wemd.rc.get_work_manager()
mode = work_manager.startup()

if work_manager.mode == work_manager.MODE_MASTER: 
    try:
        gen_istates = wemd.rc.config.get_bool('system.gen_istates', False)
        
        # Process target states
        target_states = []
        if args.tstates_from:
            target_states.extend(TargetState.states_from_file(args.tstates_from, system.pcoord_dtype))
        if args.tstates:
            tstates_strio = cStringIO.StringIO('\n'.join(args.tstates).replace(',', ' '))
            target_states.extend(TargetState.states_from_file(tstates_strio, system.pcoord_dtype))
            del tstates_strio
        
        # Process basis states
        basis_states = []
        if args.bstates_from:
            basis_states.extend(BasisState.states_from_file(args.bstates_from))
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
        
        # Check that the total probability of basis states adds to one
        tprob = sum(bstate.probability for bstate in basis_states)
        if abs(1.0 - tprob) > len(basis_states) * EPS:
            pscale = 1/tprob
            log.warning('Basis state probabilities do not add to unity; rescaling by {:g}'.format(pscale))
            for bstate in basis_states:
                bstate.probability *= pscale
        
        # Calculate progress coordinates
        sys.stdout.write('Calculating progress coordinate values for basis states.\n')
        futures = [work_manager.submit(ops.get_pcoord, propagator, basis_state) for basis_state in basis_states]
        fmap = {future: i for i,future in enumerate(futures)}
        for future in work_manager.as_completed(futures):
            updated_basis_state = future.get_result()
            basis_states[fmap[future]].pcoord = updated_basis_state.pcoord
        
        if os.path.exists(h5file):
            if args.force:
                sys.stdout.write('Deleting existing HDF5 file {!r}.\n'.format(h5file))
                os.unlink(h5file)
            else:
                sys.stderr.write('HDF5 file {!r} already exists; exiting.\n'.format(h5file))
                sys.exit(os.EX_OSFILE)
        
        # Prepare HDF5 file
        sys.stdout.write('Creating HDF5 file {!r}.\n'.format(h5file))
        data_manager.prepare_backing()        
        data_manager.save_basis_states(basis_states)
        data_manager.save_target_states(target_states)
        
        sys.stdout.write('{:d} basis state(s) present'.format(len(basis_states)))
        if wemd.rc.verbose_mode:
            sys.stdout.write(':\n')
            sys.stdout.write('{:6s}    {:12s}    {:20s}    {:20s}    {}\n'
                             .format('ID', 'Label', 'Probability', 'Aux Reference', 'Progress Coordinate'))
            for basis_state in basis_states:
                sys.stdout.write('{:<6d}    {:12s}    {:<20.14g}    {:20s}    {}\n'.
                                 format(basis_state.state_id, basis_state.label, basis_state.probability, basis_state.auxref or '',
                                        ', '.join(map(str,basis_state.pcoord))))
        sys.stdout.write('\n')
        
        sys.stdout.write('{:d} target state(s) present'.format(len(target_states)))    
        if wemd.rc.verbose_mode and target_states:
            sys.stdout.write(':\n')
            sys.stdout.write('{:6s}    {:12s}    {}\n'.format('ID', 'Label', 'Progress Coordinate'))
            for target_state in target_states:
                sys.stdout.write('{:<6d}    {:12s}    {}\n'
                                 .format(target_state.state_id, target_state.label, ','.join(map(str,target_state.pcoord))))
        sys.stdout.write('\n')
        
        # Prepare simulation
        
        
        sys.stdout.write('Preparing initial states\n')
        segs_per_state = args.segs_per_state
        segments = []
        initial_states = []
        if gen_istates:
            istate_type = InitialState.ISTATE_TYPE_GENERATED
        else:
            istate_type = InitialState.ISTATE_TYPE_BASIS
            
        for basis_state in basis_states:
            for iseg in xrange(segs_per_state):
                initial_state = data_manager.create_initial_states(0,1)[0]
                initial_state.basis_state_id =  basis_state.state_id
                initial_state.basis_state = basis_state
                initial_state.istate_type = istate_type
                segment = Segment(weight=basis_state.probability/segs_per_state,pcoord=system.new_pcoord_array(),
                                  parent_id=-(initial_state.state_id+1), wtg_parent_ids=(-(initial_state.state_id+1),),
                                  )
                initial_states.append(initial_state)
                segments.append(segment)
                
        if gen_istates:
            futures = [work_manager.submit(ops.gen_istate, propagator, initial_state.basis_state, initial_state)
                       for initial_state in initial_states]
            for future in work_manager.as_completed(futures):
                rbstate, ristate = future.get_result()
                initial_states[ristate.state_id].pcoord = ristate.pcoord
                segments[ristate.state_id].pcoord[0] = ristate.pcoord            
        else:
            for segment, initial_state in izip(segments, initial_states):
                initial_state.pcoord = basis_state.pcoord
                segment.pcoord[0] = basis_state.pcoord
                    
                
        data_manager.update_initial_states(initial_states, n_iter=0)
        
        # Make sure we have a norm of 1
        tprob = sum(segment.weight for segment in segments)
        if abs(1.0 - tprob) > len(segments) * EPS:
            pscale = 1.0/tprob
            log.warning('Weights of initial segments do not sum to unity; scaling by {:g}'.format(pscale))
            for segment in segments:
                segment.weight *= pscale
        
        
        if args.run_we:
            # TODO: what to do if something winds up in a recycling region?
            # At the moment, fail with an exception
            region_set = system.new_region_set()
            new_region_set = we_driver.run_we(region_set, segments)
            all_bins = new_region_set.get_all_bins()
        else:
            new_region_set = system.new_region_set()
            new_region_set.assign_to_bins(segments, key=Segment.initial_pcoord)
        
        all_bins = new_region_set.get_all_bins()
        bin_occupancies = numpy.array(map(operator.attrgetter('count'), all_bins))
        target_occupancies = numpy.array(map(operator.attrgetter('target_count'), all_bins))
        segments = list(new_region_set.particles)
                    
        for segment in segments:
            segment.n_iter = 1
            segment.status = Segment.SEG_STATUS_PREPARED
            assert segment.parent_id < 0
            initial_states[segment.initial_state_id].iter_used = 1
                    
        data_manager.prepare_iteration(1, segments)
        data_manager.update_initial_states(initial_states, n_iter=0)
                    
        if wemd.rc.verbose_mode:
            sys.stdout.write('\nSegments generated:\n')
            for segment in segments:
                sys.stdout.write('{!r}\n'.format(segment))
        
        
        sys.stdout.write('''
        Total bins:            {total_bins:d}
        Initial replicas:      {init_replicas:d} in {occ_bins:d} bins, total weight = {weight:g}
        Total target replicas: {total_replicas:d}
        '''.format(total_bins=len(all_bins), init_replicas=sum(bin_occupancies), occ_bins=len(bin_occupancies[bin_occupancies > 0]),
                   weight = sum(segment.weight for segment in segments), total_replicas = sum(target_occupancies)))
        
        # Send the segments over to the data manager to commit to disk            
        #data_manager.prepare_iteration(1, segments)
        data_manager.flush_backing()    
        sys.stdout.write('Simulation prepared.\n')
        work_manager.shutdown(0)
    except:
        work_manager.shutdown(4)
        raise
    
    