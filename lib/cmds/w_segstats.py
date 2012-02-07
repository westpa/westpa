from __future__ import division, print_function
from itertools import chain, starmap, izip
import os, sys, operator, argparse
import numpy
import wemd
from math import ceil, log10

import logging
log = logging.getLogger('w_segstats')

parser = wemd.rc.common_arg_parser('w_segstats', 
                                   description='''Display statistics about segments''')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
parser.add_argument('h5file', nargs='?',
                    help='WEMD HDF5 file to examine (default: determine from wemd.cfg)')
args = parser.parse_args()
output_file = args.output_file
wemd.rc.config_logging(args, tool_logger_name = 'w_segstats')

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
# column 1:    number of segments
# column 2:    continuations
# column 3:    merges
# column 4:    recycles
# column 5:    average wallclock time (s)
# column 6:    standard deviation of wallclock time (s)
# column 7:    average cpu time (s)
# column 8:    standard deviation of cpu time (s)
'''.format())
    
max_iter = sim_manager.data_manager.current_iteration - 1
for n_iter in xrange(1, max_iter+1):
    seg_index_raw = sim_manager.data_manager.get_iter_group(n_iter)['seg_index'][...]
    seg_index = numpy.empty((len(seg_index_raw),), dtype=wemd.data_manager.seg_index_dtype)
    seg_index[:] = seg_index_raw[:]
    
    n_segs = len(seg_index)
    n_continues = (seg_index['endpoint_type'] == wemd.Segment.SEG_ENDPOINT_CONTINUES).sum()
    n_merges    = (seg_index['endpoint_type'] == wemd.Segment.SEG_ENDPOINT_MERGED).sum()
    n_recycles  = (seg_index['endpoint_type'] == wemd.Segment.SEG_ENDPOINT_RECYCLED).sum()
    avg_cpu     = seg_index['cputime'].mean()
    std_cpu     = seg_index['cputime'].std()
    avg_wall    = seg_index['walltime'].mean()
    std_wall    = seg_index['walltime'].std()
    
    output_file.write(('{n_iter:10d}    {n_segs:10d}    {n_continues:10d}    {n_merges:10d}    {n_recycles:10d}    '
                      +'{avg_wall:12.4g}    {std_wall:12.4g}    {avg_cpu:12.4g}    {std_cpu:12.4g}\n')
                      .format(**locals()))
