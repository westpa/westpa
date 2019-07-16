import os
import shutil
import argparse


def remove_traj_segs(basedir, first_iter, last_iter, iter_prec):
    '''Delete the trajectory segment directions between first_iter and
    last_iter. Looks for directories in basedir'''

    iter_template = '{:0' + str(iter_prec) + 'd}'

    for k in range(first_iter, last_iter):
        dirname = os.path.join(basedir, iter_template.format(k))
        if os.path.exists(dirname):
            print(dirname)
            shutil.rmtree(dirname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--first_iter', required=True, type=int, help='First iteration to remove')
    parser.add_argument('--last_iter', required=True, type=int, help='Final iteration to remove')
    parser.add_argument('-D', '--dir', default='traj_segs', help='Base directory containing trajectory seg directories')
    parser.add_argument('-i', '--iter_prec', default=6, type=int, help='Width of iteration for padding')

    args = parser.parse_args()

    remove_traj_segs(args.dir, args.first_iter, args.last_iter, args.iter_prec)
