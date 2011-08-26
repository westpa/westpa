from __future__ import division, print_function

import os, sys, logging, numpy, operator, argparse
log = logging.getLogger('w_init')

import wemd

parser = argparse.ArgumentParser('w_init', description='initialize a new WEMD simulation')
wemd.rc.add_common_args(parser)
parser.add_argument('--force', dest='force', action='store_true',
                         help='overwrite any existing simulation data')
parser.add_argument('--ptol', dest='ptol', type=float, default=1.0e-8,
                         help='tolerance for sum of initial/recycle probabilities (default: 1.0e-8)')

(args, aux_args) = parser.parse_known_args()

if aux_args:
    sys.stderr.write('unexpected command line argument(s) encountered: {}\n'.format(aux_args))
    sys.exit(os.EX_USAGE)

wemd.rc.process_common_args(args)
system = wemd.rc.get_system_driver()
h5file = wemd.rc.config.get_path('data.h5file')
data_manager = wemd.rc.get_data_manager()
data_manager.backing_file = h5file
data_manager.system = system

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

# Prepare simulation
system.prepare_run()
region_set = system.region_set

tiprob = 0.0
trprob = 0.0
sys.stdout.write('\nInitial state:\n')
sys.stdout.write('{:<16} {:<12} {:<12} {}\n'.format('Label', 'Init. Prob.', 'Recyc. Prob.', 'Coordinates'))
for istate in system.initial_states:
    sys.stdout.write('{istate.label:<16} {istate.initial_prob:<12g} {istate.recycle_prob:<12g} {pcoord!s:<52}\n'\
                     .format(istate=istate, pcoord=list(istate.pcoord)))
    tiprob += istate.initial_prob
    trprob += istate.recycle_prob

MACHEPS = numpy.finfo(numpy.float64).eps
if abs(1.0 - tiprob) > args.ptol:
    sys.stderr.write('Initial probabilities do not sum to one.')
    sys.exit(1)
if abs(1.0 - trprob) > args.ptol:
    sys.stderr.write('Recycle probabilities do not sum to one.')
    sys.exit(1)

# Create initial segments
segments = []
for (i_istate, istate) in enumerate(system.initial_states):
    # Skip microstates that are for recycling only
    if istate.initial_prob == 0.0: continue
    target_count = istate.bin.target_count
    for i in xrange(0, target_count):
        segment = wemd.Segment(pcoord = system.new_pcoord_array(),
                               weight = istate.initial_prob / target_count,
                               p_parent_id = -(i_istate+1),
                               parent_ids = set([-(i_istate+1)]),
                               status = wemd.Segment.SEG_STATUS_PREPARED)
        segment.pcoord[0] = istate.pcoord
        istate.bin.add(segment)
        segments.append(segment)
    sys.stdout.write('%d replicas from initial point %r\n' % (target_count,istate.label))

iprobtot = region_set.weight
all_bins = region_set.get_all_bins()
bin_occupancies = numpy.array(map(operator.attrgetter('count'), all_bins))
target_occupancies = numpy.array(map(operator.attrgetter('target_count'), all_bins))

sys.stdout.write('''
Total bins:            {total_bins:d}
Initial replicas:      {init_replicas:d} in {occ_bins:d} bins, total weight = {weight:g}
Total target replicas: {total_replicas:d}
'''.format(total_bins=len(all_bins), init_replicas=sum(bin_occupancies), occ_bins=len(bin_occupancies[bin_occupancies > 0]),
           weight = iprobtot, total_replicas = sum(target_occupancies)))

# The user-side check for this was above; this is an assertion that the above assignment to bins 
# and division of probability is correct
assert abs(system.region_set.weight - tiprob) < MACHEPS*sum(bin_occupancies)

# Send the segments over to the data manager to commit to disk            
data_manager.prepare_iteration(1, segments)
data_manager.flush_backing()

system.region_set.clear()    
sys.stdout.write('Simulation prepared.\n')
