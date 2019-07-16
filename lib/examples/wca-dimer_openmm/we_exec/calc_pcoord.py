import numpy as np
import tables
import argparse
import os
import sys


def run_npz(file, out):
    data = np.load(file)
    coords = data['coord']

    x = coords[0,:]
    y = coords[1,:]

    d = 10.0 * np.sqrt(np.sum((x - y)**2))
    print(d)

    with open(out, 'w') as f:
        f.write('{:f}\n'.format(d))


def run_h5(file, out):
    f = tables.File(file, 'r')
    coords = f.root.coordinates[:]
    print('coords shape: ', coords.shape)

    x = coords[:,0,:]
    y = coords[:,1,:]

    d = 10.0*np.sqrt(np.sum((x - y)**2, axis=1))

    np.savetxt(out, d)

    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='''Name of input file (either
                a mdtraj HDFReporter file or a .npz file''')
    parser.add_argument('-o', '--output', help='output file')

    args = parser.parse_args()

    if args.file.endswith('.h5'):
        run_h5(args.file, args.output)
    elif args.file.endswith('.npz'):
        run_npz(args.file, args.output)
    else:
        sys.exit(1)
