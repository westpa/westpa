from __future__ import division, print_function
import os, sys
from math import ceil, log10

import logging
log = logging.getLogger('w_binprobs')

import wemd

parser = wemd.rc.common_arg_parser()
parser.add_argument('-o', '--output', dest='output_file',
                    help='Store output in OUTPUT_FILE (default: write to standard output).')
parser.add_argument('-l', '--labels', dest='print_labels', action='store_true',
                    help='Print bin labels (to stdout) corresponding to each column in the output '
                        +'(requires system driver to be available; default: do not print labels)')
parser.add_argument('-p', '--precision', dest='precision', type=int, 
                    help='Numer of significant figures for probability display (default: 6)',
                    default=6)
args = parser.parse_args()

wemd.rc.config_logging(args)
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update({'args.%s' % key : value for (key,value) in args.__dict__.viewitems() if not key.startswith('_')})
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()

if args.print_labels:
    sim_manager.load_system_driver()
    bins = sim_manager.system.region_set.get_all_bins()
    maxwidth = int(ceil(log10(len(bins))))
    sys.stdout.write('bin labels:\n')
    for (ibin, bin) in enumerate(bins):
        sys.stdout.write('  column {0:<{maxwidth}d}: {1!s}\n'.format(ibin, bin.label, maxwidth=maxwidth))

if args.output_file:
    output_file = open(args.output_file, 'wt')
else:
    output_file = sys.stdout
    sys.stdout.write('----\n')

max_iter = sim_manager.data_manager.current_iteration-1
max_iter_width = int(ceil(log10(max_iter+1)))
prob_width = min(12,args.precision + 5)
for iiter in xrange(0, max_iter):
    iter_group = sim_manager.data_manager.get_iter_group(iiter+1)
    prob_fields = ['{:<{width}.{precision}g}'.format(prob, width=prob_width, precision=args.precision) 
                   for prob in iter_group['bin_probs'][...]]
    output_file.write('{0:<{maxwidth}d} {1}\n'.format(iiter+1, ' '.join(prob_fields), maxwidth = max_iter_width))
