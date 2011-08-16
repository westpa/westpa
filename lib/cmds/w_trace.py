from __future__ import division, print_function
import sys, argparse
import numpy
import wemd, wemdtools

import logging
log = logging.getLogger('w_trace')

parser = wemd.rc.common_arg_parser('w_trace', 
                                   description='''Trace trajectories''')
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).',
                    type=argparse.FileType('wt'), default=sys.stdout)
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
parser.add_argument('-i', '--input', '--datafile', dest='datafile',
                    help='WEMD HDF5 file to examine (default: determine from wemd.cfg)')
parser.add_argument('segment', nargs='+',
                    help='Segments to trace, specified as "n_iter:seg_id"')
args = parser.parse_args()
output_file = args.output_file
wemd.rc.config_logging(args, tool_logger_name = 'w_trace')

runtime_config = None
sim_manager = None
data_manager = None
if args.datafile:
    # Reading some other HDF5 file
    data_manager = wemd.data_manager.WEMDDataManager(sim_manager = None, backing_file = args.datafile)
else:
    # Reading from HDF5 file specified in wemd.cfg
    runtime_config = wemd.rc.read_config(args.run_config_file)
    runtime_config.update_from_object(args)
    sim_manager = wemd.rc.load_sim_manager(runtime_config)
    sim_manager.load_data_manager()
    data_manager = sim_manager.data_manager

data_manager.open_backing(mode='r')
if not args.suppress_headers:
    output_file.write('''\
# column 0:    terminal segment (formatted as "n_iter:seg_id")
# column 1:    n_iter
# column 2:    seg_id 
''')

ttree = wemdtools.trajectories.trajtree.TrajTree(data_manager)
for segspec in args.segment:
    n_iter, seg_id = map(int, segspec.split(':'))
    trajectory = ttree.trace_to_root(n_iter, seg_id)
    osegspec = '{n_iter:d}:{seg_id:d}'.format(n_iter = n_iter, seg_id = seg_id)
    for segment in trajectory:
        output_file.write('{osegspec:20s}    {segment.n_iter:10d}    {segment.seg_id:10d}\n'
                          .format(osegspec = osegspec, segment = segment))     
