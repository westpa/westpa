from __future__ import division, print_function
import os, sys

import logging
log = logging.getLogger('w_binprobs')

import wemd

parser = wemd.rc.common_arg_parser()
args = parser.parse_args()

wemd.rc.config_logging(args)
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update({'args.%s' % key : value for (key,value) in args.__dict__.viewitems() if not key.startswith('_')})
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
sim_manager.data_manager.open_backing()

for iiter in xrange(0, sim_manager.data_manager.current_iteration-1):
    iter_group = sim_manager.data_manager.get_iter_group(iiter+1)
    prob_fields = ['{:20.16g}'.format(prob) for prob in iter_group['bin_probs'][...]]
    sys.stdout.write('{0:8d} {1}\n'.format(iiter+1, ' '.join(prob_fields)))
