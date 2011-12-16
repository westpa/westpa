from __future__ import division, print_function
import sys, argparse
import numpy
import wemd

import logging
log = logging.getLogger('w_succ')

parser = wemd.rc.common_arg_parser('w_succ', 
                                   description='''List segments which successfully reach a target state''')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
parser.add_argument('h5file', nargs='?',
                    help='WEMD HDF5 file to examine (default: determine from wemd.cfg)')
args = parser.parse_args()
output_file = args.output_file
wemd.rc.config_logging(args, tool_logger_name = 'w_succ')

if args.h5file:
    from wemd.util.config_dict import ConfigDict
    runtime_config = ConfigDict()
    runtime_config['data.h5file'] = args.h5file
else:
    # Load from wemd.cfg
    runtime_config = wemd.rc.read_config(args.run_config_file)
    
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()

if not args.suppress_headers:
    output_file.write('''\
# successful (recycled) segments 
# column 0:    iteration
# column 1:    seg_id
# column 2:    weight
'''.format())
    
max_iter = sim_manager.data_manager.current_iteration - 1
for n_iter in xrange(1, max_iter+1):
    seg_index_raw = sim_manager.data_manager.get_iter_group(n_iter)['seg_index'][...]
    seg_index = numpy.empty((len(seg_index_raw),), dtype=wemd.data_manager.seg_index_dtype)
    seg_index[:] = seg_index_raw[:]
    for (seg_id, seg_info) in enumerate(seg_index):
        if seg_info['endpoint_type'] == wemd.Segment.SEG_ENDPOINT_RECYCLED:
            output_file.write('{n_iter:10d}    {seg_id:10d}    {weight!r:<24s}\n'
                              .format(n_iter = n_iter, seg_id = seg_id, weight = seg_info['weight']))
