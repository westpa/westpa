from __future__ import division, print_function
from itertools import chain, starmap, izip
import os, sys, operator, argparse
import numpy
import wemd
from math import ceil, log10

import logging
log = logging.getLogger('w_dumpsegs')

parser = wemd.rc.common_arg_parser('w_dumpsegs', description='''\
Dump segment data for a given iteration. Useful for debugging the data manager.''')
parser.add_argument('-p', '--print-pcoords', dest='print_pcoords', action='store_true',
                    help='print initial and final progress coordinates for each segment')
parser.add_argument('-i', '--iteration', dest='n_iter', type=int,
                    help='Use data from iteration N_ITER (default: last complete iteration)')
parser.add_argument('-s', '--segments', dest='seg_ids', type=long, nargs='+',
                    help='Print only the given SEG_IDS')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
args = parser.parse_args()
output_file = args.output_file

wemd.rc.config_logging(args)
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()
#sim_manager.load_system_driver()
#bins = sim_manager.system.region_set.get_all_bins()

n_iter = args.n_iter or sim_manager.data_manager.current_iteration - 1
segments = sim_manager.data_manager.get_segments(n_iter)
seg_id_pool = set(len(str(segment.seg_id)) for segment in segments)
seg_id_pool.update(len(str(segment.p_parent_id)) for segment in segments)

max_seg_id_len = max(seg_id_pool)
max_status_name_len = max(map(len,wemd.Segment.status_names.values()))
max_endpoint_type_len = max(map(len, wemd.Segment.endpoint_type_names.values()))
max_n_parents_len = len(str(max(segment.n_parents for segment in segments)))

report_line = ( '{segment.n_iter:d}  {segment.seg_id:{max_seg_id_len}d}  {segment.weight:20.14g}' 
                +'  {status_name:{max_status_name_len}s} ({segment.status})'
                +'  {segment.walltime:<12.6g} {segment.cputime:<12.6g}'
                +'  {endpoint_type_name:{max_endpoint_type_len}s} ({segment.endpoint_type})'
                +'  {segment.n_parents:{max_n_parents_len}d} {segment.p_parent_id:{max_seg_id_len}d} {parents_str}' 
                +'\n')
pcoord_lines = ('  pcoord[0]  = {init_pcoord}\n  pcoord[-1] = {final_pcoord}'
                +'\n')

seg_ids = set(args.seg_ids or [])
for (seg_id, segment) in enumerate(segments):
    if seg_ids and seg_id not in seg_ids: continue
    parents_str = '['+', '.join(map(str,sorted(segment.parent_ids)))+']'
    init_pcoord_str = '[' + ', '.join('{pcval:<12.6g}'.format(pcval=pce) for pce in segment.pcoord[0]) + ']'
    final_pcoord_str = '[' + ', '.join('{pcval:<12.6g}'.format(pcval=pce) for pce in segment.pcoord[-1]) + ']'
    output_file.write(report_line.format(segment=segment, 
                                         status_name = segment.status_names[segment.status],
                                         endpoint_type_name = segment.endpoint_type_names[segment.endpoint_type],
                                         parents_str = parents_str,
                                         max_seg_id_len=max_seg_id_len,
                                         max_status_name_len=max_status_name_len,
                                         max_endpoint_type_len=max_endpoint_type_len,
                                         max_n_parents_len=max_n_parents_len))
    if args.print_pcoords:
        output_file.write(pcoord_lines.format(init_pcoord = init_pcoord_str,final_pcoord = final_pcoord_str))
