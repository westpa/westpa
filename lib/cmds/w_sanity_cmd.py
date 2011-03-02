from __future__ import division, print_function
from itertools import chain, starmap, izip
import os, sys, operator, argparse
import numpy
import wemd
from math import ceil, log10

import logging
log = logging.getLogger('w_sanity')

parser = wemd.rc.common_arg_parser('w_sanity', description='''Check a simulation data file for consistency.''')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
args = parser.parse_args()
output_file = args.output_file

wemd.rc.config_logging(args, tool_logger_name = 'w_sanity')
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)

sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()
#sim_manager.load_system_driver()
#bins = sim_manager.system.region_set.get_all_bins()


prev_segs = []
curr_segs  = []
next_segs   = []

n_errors = 0
for n_iter in xrange(1, sim_manager.data_manager.current_iteration):
    output_file.write('examining iteration {}\n'.format(n_iter))
    if n_iter > 1:
        prev_segs = curr_segs
        curr_segs = next_segs
    else:
        curr_segs = sim_manager.data_manager.get_segments(1)
        
    if n_iter < sim_manager.data_manager.current_iteration:
        next_segs = sim_manager.data_manager.get_segments(n_iter+1)
        
    prev_segs_by_id = {segment.seg_id: segment for segment in prev_segs}
    curr_segs_by_id = {segment.seg_id: segment for segment in curr_segs}
    next_segs_by_id = {segment.seg_id: segment for segment in next_segs}
        
    # Check to see if seg_ids are correct
    for (iseg, segment) in enumerate(curr_segs):
        if segment.seg_id != iseg:
            n_errors += 1
            output_file.write('segment {seg_id:d}: seg_id does not match location in HDF5 file ({iseg:d})\n'
                              .format(seg_id=seg_id, iseg=iseg))
    
    # Check to see if current iteration's start positions are the same as the previous iteration's end position
    if prev_segs:
        for (iseg, segment) in enumerate(curr_segs):
            if segment.status == segment.SEG_STATUS_COMPLETE and segment.p_parent_id > 0:
                try:
                    parent = prev_segs_by_id[segment.p_parent_id]
                except KeyError:
                    n_errors += 1
                    output_file.write('segment {segment.seg_id} has invalid p_parent_id ({segment.p_parent_id})\n'
                                      .format(segment=segment))
                else:
                    assert parent.n_iter == n_iter - 1
                    parent_final_pcoord = parent.pcoord[-1]
                    current_initial_pcoord = segment.pcoord[0]
                    if not (current_initial_pcoord == parent_final_pcoord).all():
                        n_errors += 1
                        output_file.write(('initial pcoord {init_pcoord!s} of segment {segment.seg_id}'
                                           +' does not match final pcoord {final_pcoord!s} of parent {segment.p_parent_id}\n')
                                           .format(segment=segment, parent=parent,
                                                   init_pcoord=current_initial_pcoord,
                                                   final_pcoord=parent_final_pcoord))
                     
output_file.write('{} errors found\n'.format(n_errors))
if n_errors:
    sys.exit(1)
else:
    sys.exit(0)
    
    
    

