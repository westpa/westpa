from __future__ import division, print_function
from itertools import chain, starmap, izip
import os, sys, operator, argparse
import numpy
import wemd
from math import ceil, log10

import logging
log = logging.getLogger('w_ntop')

parser = wemd.rc.common_arg_parser('w_ntop', description='''Retrieve a number of high-weight replicas from each bin.''')
parser.add_argument('-N', '--perbin', dest='perbin', type=int,
                    help='Retrieve the PERBIN highest weighted replicas from each bin (default: 1)',
                    default=1)
parser.add_argument('-i', '--iteration', dest='n_iter', type=int,
                    help='Use data from iteration N_ITER (default: last complete iteration)')
parser.add_argument('-l', '--labels', dest='print_labels', action='store_true',
                    help='Print bin labels (to stdout) corresponding to each column in the output '
                        +'(default: do not print labels)')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
parser.add_argument('--output-istates', dest='istates_file',
                    help='Write an initial states file to ISTATES_FILE (default: do not write).',
                    type=argparse.FileType('wt'), default=None)
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
args = parser.parse_args()
output_file = args.output_file
istates_file = args.istates_file


wemd.rc.config_logging(args)
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()
sim_manager.load_system_driver()
bins = sim_manager.system.region_set.get_all_bins()

if args.print_labels:
    maxwidth = int(ceil(log10(len(bins))))
    sys.stdout.write('bin labels:\n')
    for (ibin, bin) in enumerate(bins):
        sys.stdout.write('  bin {0:<{maxwidth}d}: {1!s}\n'.format(ibin, bin.label, maxwidth=maxwidth))


n_replicas = args.perbin
max_iter = sim_manager.data_manager.current_iteration
n_iter = args.n_iter or max_iter
if n_iter > max_iter:
    sys.stderr.write('invalid iteration specified (max iteration is {:d}\n'.format(max_iter))
    sys.exit(1)

sys.stdout.write('reporting on iteration {:d}\n'.format(n_iter))
sys.stdout.write('will report the {0:d} highest-probability replica(s) for each of {1:d} bins\n'.format(n_replicas,
                                                                                                        len(bins)))
if output_file is sys.stdout:
    sys.stdout.write('----\n')

# Get segments
segments = sim_manager.data_manager.get_segments(n_iter)

# Bin the starting positions of each segment for the given iteration
pcoords = numpy.empty((len(segments),) + segments[0].pcoord[0].shape, segments[0].pcoord.dtype)
for (iseg, segment) in enumerate(segments):
    pcoords[iseg] = segment.pcoord[0]
bin_assignments = sim_manager.system.region_set.map_to_bins(pcoords)
assert len(bin_assignments) == len(segments)
for (iseg, bin) in enumerate(bin_assignments):
    bin.add(segments[iseg])

# bin index, index_in_bin, seg_id, weight, pcoord
ntop = []
for (ibin, bin) in enumerate(bins):
    sorted_segs = list(sorted(bin, key=operator.attrgetter('weight')))
    for (iseg, segment) in enumerate(sorted_segs[:n_replicas]):
        ntop.append( (ibin, iseg, long(segment.seg_id), long(segment.p_parent_id), segment.weight, tuple(segment.pcoord[0]) ) )

max_ibin_len = int(ceil(log10(len(bins))))
#max_seg_id_len = int(ceil(log10(abs(max(segment.seg_id for segment in segments)))))
max_seg_id_len = int(ceil(log10(max(max(segment.seg_id, abs(segment.p_parent_id or 0)) for segment in segments)))) + 1
max_itop_len = int(ceil(log10(n_replicas)))
total_prob = sum(top[4] for top in ntop)

if not args.suppress_headers:
    output_file.write('''\
# {n_replicas:d} highest-weight particle(s) for each bin 
# at beginning of iteration {n_iter:d}
# column 0:    bin index
# column 1:    weight rank of segment within bin
# column 2:    seg_id
# column 3:    p_parent_id
# column 4:    weight
# columns >=5: pcoord
'''.format(n_replicas = n_replicas, n_iter = n_iter))
    
    if istates_file is not None:
        istates_file.write('''\
# column 0:    state label
# column 1:    initial probability
# column 2:    recycling probability
# columns >=3: pcoord
''')

for top in ntop:
    (ibin, itop, seg_id, p_parent_id, weight, pcoord) = top
    pcoord_txt = ' '.join('{!r:<24s}'.format(q) for q in pcoord)
    output_file.write(( '{ibin:{max_ibin_len}d}  {itop:{max_itop_len}d}  {seg_id:<{max_seg_id_len}d}'
                       +'  {p_parent_id:<{max_seg_id_len}d}  {weight!r:<24s}  {pcoord_txt}\n')\
                      .format(ibin=ibin, itop=itop+1, seg_id=seg_id, weight=weight, pcoord_txt=pcoord_txt,
                              p_parent_id=p_parent_id,
                              max_ibin_len=max_ibin_len, max_seg_id_len=max_seg_id_len, max_itop_len=max_itop_len))

    # If requested, write an initial states file
    if istates_file is not None:
        rel_prob = weight/total_prob
        state_label = 'bin_{:d}'.format(ibin)
        istates_file.write('{state_label:<24s}  {prob!r:<24s}  {prob!r:<24s}  {pcoord_txt}\n'\
                           .format(state_label=state_label, prob=rel_prob, pcoord_txt=pcoord_txt ))

