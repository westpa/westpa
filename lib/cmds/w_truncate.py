from __future__ import division, print_function

import argparse

import logging
log = logging.getLogger('w_truncate')

import westpa

parser = argparse.ArgumentParser('w_truncate', description='''\
Remove all iterations after a certain point in a  
''')
westpa.rc.add_args(parser)
parser.add_argument('-n', '--iter', dest='n_iter', type=int,
                    help='Truncate this iteration and those following.')
args = parser.parse_args()
westpa.rc.process_args(args, config_required=False)
dm = westpa.rc.get_data_manager()
dm.open_backing()
max_iter = dm.current_iteration
n_iter = args.n_iter if args.n_iter > 0 else dm.current_iteration

for i in xrange(n_iter, dm.current_iteration+1):
    dm.del_iter_group(i)

dm.del_iter_summary(n_iter)
dm.current_iteration = n_iter - 1

print('simulation data truncated after iteration {}'.format(dm.current_iteration))

dm.flush_backing()
dm.close_backing()