import logging
from argparse import ArgumentParser

import westpa


log = logging.getLogger('w_truncate')

warning_string = '''\
NOTE: w_truncate only deletes iteration groups from the HDF5 data store.
It is recommended that any iteration data saved to the file system (e.g. in the
traj_segs directory) is deleted or moved for the corresponding iterations.
'''


def entry_point():
    parser = ArgumentParser(
        'w_truncate',
        description='''\
    Remove all iterations after a certain point in a WESTPA simulation.
    ''',
        epilog=warning_string,
    )

    westpa.rc.add_args(parser)
    parser.add_argument(
        '-W',
        '--west-data',
        dest='we_h5filename',
        metavar='WEST_H5FILE',
        help='''Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in west.cfg).''',
    )
    parser.add_argument('-n', '--iter', dest='n_iter', type=int, default=0, help='Truncate this iteration and those following.')

    args = parser.parse_args()
    westpa.rc.process_args(args, config_required=False)
    dm = westpa.rc.get_data_manager()
    if args.we_h5filename:
        dm.we_h5filename = args.we_h5filename

    dm.open_backing()
    # max_iter = dm.current_iteration

    if args.n_iter > dm.current_iteration:
        parser.error(
            'Provided iteration {} > current iteration {} of the {} HDF5 file. Exiting without doing anything.'.format(
                args.n_iter, dm.current_iteration, dm.we_h5filename.split('/')[-1]
            )
        )
    else:
        n_iter = args.n_iter if args.n_iter > 0 else dm.current_iteration

    for i in range(n_iter, dm.current_iteration + 1):
        dm.del_iter_group(i)

    dm.del_iter_summary(n_iter)
    dm.current_iteration = n_iter - 1

    westpa.rc.pstatus('simulation data truncated after iteration {}'.format(dm.current_iteration))
    westpa.rc.pstatus('\n' + warning_string)
    westpa.rc.pflush()

    dm.flush_backing()
    dm.close_backing()


if __name__ == '__main__':
    entry_point()
