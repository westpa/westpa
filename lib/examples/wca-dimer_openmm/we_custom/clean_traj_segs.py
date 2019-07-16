import os
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--iter-start', default=1, type=int)
parser.add_argument('--iter-stop', required=True, type=int)
parser.add_argument('-n', '--dry-run', action='store_true', default=False)

args = parser.parse_args()

for k in range(args.iter_start, args.iter_stop + 1):
    files = glob.glob('traj_segs/iter_{:06d}_*.npz'.format(k))
    if args.dry_run:
        print(files)
    else:
        for f in files:
            if not args.dry_run:
                os.remove(f)

