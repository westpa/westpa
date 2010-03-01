import os, sys, re, tempfile
from itertools import chain
import numpy
import h5py
from optparse import OptionParser

parser = OptionParser(usage='%prog [OPTIONS] TEXT_FILE HDF5_FILE HDF5_NODE',
                      description = 'Import data from a text file into an HDF5 file')
parser.add_option('-t', '--timestep', dest='dt',
                  help='simulation timestep (default: read from input)')
parser.add_option('-d', '--ndim', dest='ndim', type='int',
                  help='number of columns '
                      +'(default: read from input)')
(opts, args) = parser.parse_args()

if len(args) != 3:
    sys.stderr.write('exactly 3 arguments are required\n')
    parser.print_help(sys.stderr)
    sys.exit(2)
    
input_filename = args[0]
hdf5_filename = args[1]
hdf5_nodename = args[2]

if input_filename == '-':
    input_file = sys.stdin
else:
    input_file = open(input_filename, 'rt')

h5file = h5py.File(hdf5_filename)

dt = opts.dt
ndim = opts.ndim
linebuffer = []
if dt is None or ndim is None:
    line1 = input_file.readline()
    while line1 and line1[0] in ('#', '@'):
        line1 = input_file.readline()
        continue
    linebuffer.append(line1)
    line2 = input_file.readline()
    linebuffer.append(line2)
    fields1 = line1.split()
    fields2 = line2.split()
    if dt is None:
        dt = float(fields2[0]) - float(fields1[0])
        sys.stdout.write('using dt = %g\n' % dt)
    if ndim is None:
        ndim = len(fields1)-1
        sys.stdout.write('using %d-dimensional coordinate\n' % ndim)

chunk_size = 32768
conv_buffer = numpy.empty((chunk_size,ndim), numpy.float64)
chunks_read = 0
lines_read = 0

nodenames = hdf5_nodename.split('/')
grp = None
for nodename in nodenames[:-1]:
    grp = (grp or h5file).create_group(nodename)

dset = (grp or h5file).create_dataset(nodenames[-1], (1,ndim), numpy.float64,
                                      maxshape=(None, ndim),
                                      chunks=(chunk_size,ndim)) 
dset.attrs['timestep'] = dt
for line in chain(linebuffer, input_file):
    if line[0] in ('#', '@'): continue
    
    fields = line.split()
    lines_read += 1
    if ndim == 1:
        conv_buffer[lines_read-1,0] = float(fields[1])
    else:
        conv_buffer[lines_read-1,:] = [float(field) for field in fields[1:ndim+1]]
    if lines_read == chunk_size:
        new_shape = (chunk_size*(chunks_read+1),) + tuple(dset.shape[1:])
        dset.resize(new_shape)
        dset[-chunk_size:,0,:] = conv_buffer
        lines_read = 0
        chunks_read += 1
        
if lines_read:
    old_shape = dset.shape
    new_shape = (dset.shape[0] + lines_read,) + dset.shape[1:]
    dset.resize(new_shape)
    dset[-lines_read:,0,:] = conv_buffer[:lines_read]
        
h5file.close()
input_file.close()
    

