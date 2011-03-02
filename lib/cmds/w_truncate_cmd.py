from __future__ import division, print_function
import os, sys, argparse

import logging
log = logging.getLogger('w_truncate')

import wemd

parser = wemd.rc.common_arg_parser()
parser.add_argument('-n', '--iter', dest='iter', default=None, help='Truncate after this iteration (by default remove the last iteration)')

args = parser.parse_args()

wemd.rc.config_logging(args)
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
dm = sim_manager.data_manager
dm.open_backing()
max_iter = sim_manager.data_manager.current_iteration

if args.iter is None:
    iter = max_iter 
else:
    iter = int(args.iter)

for i in xrange(iter, max_iter + 1):
    dm.del_iter(i)

dm.del_iter_summary(iter)
sim_manager.data_manager.current_iteration = iter - 1

dm.flush_backing()
dm.close_backing()