


import argparse

import logging
log = logging.getLogger('w_truncate')

import westpa

warning_string = '''\
NOTE: w_truncate only deletes iteration groups from the HDF5 data store.
It is recommended that any iteration data saved to the file system (e.g. in the
traj_segs directory) is deleted or moved for the corresponding iterations.
'''

parser = argparse.ArgumentParser('w_truncate', description='''\
Remove all iterations after a certain point in a WESTPA simulation.
''',
epilog=warning_string)


westpa.rc.add_args(parser)
parser.add_argument('-n', '--iter', dest='n_iter', type=int,
                    help='Truncate this iteration and those following.')
args = parser.parse_args()
westpa.rc.process_args(args, config_required=False)
dm = westpa.rc.get_data_manager()
dm.open_backing()
max_iter = dm.current_iteration
n_iter = args.n_iter if args.n_iter > 0 else dm.current_iteration

for i in range(n_iter, dm.current_iteration+1):
    dm.del_iter_group(i)

dm.del_iter_summary(n_iter)
dm.current_iteration = n_iter - 1

print('simulation data truncated after iteration {}'.format(dm.current_iteration))
print('\n' + warning_string)

dm.flush_backing()
dm.close_backing()
