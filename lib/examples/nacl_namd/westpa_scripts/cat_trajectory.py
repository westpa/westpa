#!/usr/bin/python

import h5py, numpy, sys

infile = numpy.loadtxt(sys.argv[1], usecols = (0, 1))
west   = h5py.File('west.h5')
coords = []
for iteration, seg_id in infile[1:]:
    iter_key = "iter_{0:08d}".format(int(iteration))
    SOD      = west['iterations'][iter_key]['auxdata']['coord'][seg_id,1:,0,:]
    CLA      = west['iterations'][iter_key]['auxdata']['coord'][seg_id,1:,1,:]
    coords  += [numpy.column_stack((SOD, CLA))]
with open(sys.argv[1][:-4] + ".xyz", 'w') as outfile:
    for i, frame in enumerate(numpy.concatenate(coords)):
        outfile.write("2\n")
        outfile.write("{0}\n".format(i))
        outfile.write("SOD {0:9.5f} {1:9.5f} {2:9.5f}\n".format(
          float(frame[0]), float(frame[1]), float(frame[2])))
        outfile.write("CLA {0:9.5f} {1:9.5f} {2:9.5f}\n".format(
          float(frame[3]), float(frame[4]), float(frame[5])))
