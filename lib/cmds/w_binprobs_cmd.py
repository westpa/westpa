from __future__ import division, print_function
import os, sys, argparse, numpy
from math import ceil, log10

import logging
log = logging.getLogger('w_binprobs')

import wemd, wemdtools

parser = wemd.rc.common_arg_parser(description='''\
Calculate the instantaneous probability in each bin at the beginning of each iteration.
''')
# Region set options
wemdtools.bins.add_region_set_options(parser)
# Output options
parser.add_argument('-o', '--output', dest='output_file', type=argparse.FileType('wt'), default=sys.stdout,
                    help='Store per-iteration output in OUTPUT_FILE (default: write to standard output).')
parser.add_argument('-l', '--labels', dest='print_labels', action='store_true',
                    help='Print bin labels corresponding to each column in the output '
                        +'(default: do not print labels)')
parser.add_argument('-p', '--precision', dest='precision', type=int, 
                    help='Number of significant figures for probability display (default: 6)',
                    default=6)
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')

args = parser.parse_args()
wemd.rc.config_logging(args, 'w_binprobs')
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()

region_set = wemdtools.bins.get_region_set_from_args(args, status_stream=sys.stdout)
bins = region_set.get_all_bins()
bin_map = {id(bin): ibin for (ibin, bin) in enumerate(bins)}
n_bins = len(bins)

max_iter = sim_manager.data_manager.current_iteration
max_iter_width = int(ceil(log10(max_iter+1)))
prob_width = min(12,args.precision + 5)

if not args.suppress_headers:
    args.output_file.write('''\
# WE bin population analysis
# ----
''')

if args.print_labels:
    maxwidth = int(ceil(log10(n_bins)))
    args.output_file.write('# bin labels:\n')
    for (ibin, bin) in enumerate(bins):
        args.output_file.write('#  bin {0:<{maxwidth}d}: {1!s}\n'.format(ibin, bin.label, maxwidth=maxwidth))
    args.output_file.write('# ----\n')
args.output_file.write('''\
# column 0:  iteration
# following: bin population at beginning of iteration
# ----
''')

# Retrieve bin probabilities
binprobs = numpy.empty((n_bins,), numpy.float64)
for i_iter in xrange(0, max_iter):
    n_iter = 1 + i_iter
    binprobs[:] = 0.0
    try:
        iter_group = sim_manager.data_manager.get_iter_group(n_iter)
    except KeyError:
        # current_iteration is set but segments don't exist
        break
    
    weights = iter_group['seg_index'][:, 'weight']

    # Read only the initial point
    pcoords = iter_group['pcoord'][:,0,:]
    bin_indices = [bin_map[id(bin)] for bin in region_set.map_to_bins(pcoords)]
    for (seg_id, ibin) in enumerate(bin_indices):
        binprobs[ibin] += weights[seg_id]
        
    prob_fields = ['{:<{width}.{precision}g}'.format(prob, width=prob_width, precision=args.precision) 
                   for prob in binprobs]
    args.output_file.write('{0:<{maxwidth}d} {1}\n'.format(n_iter, ' '.join(prob_fields), maxwidth = max_iter_width))
