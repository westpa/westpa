import os, sys, re, tempfile
from itertools import chain
import numpy
import h5py
from optparse import OptionParser

parser = OptionParser(usage='%prog [OPTIONS] TEXT_FILE',
                      description = 'Create a pcoord file from a text file')
parser.add_option('-o', '--output', dest='output', default='dist.h5',
                  help='destination HDF5 file (default: dist.h5)')
parser.add_option('-t', '--timestep', dest='dt',
                  help='simulation timestep (default: read from input)')
parser.add_option('-d', '--ndim', dest='ndim', type='int',
                  help='number of progress coordinate dimensions '
                      +'(default: read from input)')
(opts, args) = parser.parse_args()

if len(args) != 1:
    sys.stderr.write('exactly one argument is required\n')
    parser.print_help(sys.stderr)
    sys.exit(2)

if args[0] == '-':
    input_file = sys.stdin
else:
    input_file = open(args[0], 'rt')
h5file = h5py.File(opts.output, 'w')

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
grp = h5file.create_group('pcoords')
dset = grp.create_dataset('pcoord', (1,1,ndim), numpy.float64,
                             maxshape=(None, None, ndim),
                             chunks=(chunk_size,1,ndim)) 
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
    

