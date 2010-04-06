import sys
from itertools import chain
import numpy
import h5py
from optparse import OptionParser

parser = OptionParser(usage='%prog [OPTIONS] TEXT_FILE HDF5_FILE HDF5_NODE',
                      description = 'Import data from a text file into an HDF5 file')
parser.add_option('--ncol', dest='ncol', type='int',
                  help='number of columns '
                      +'(default: read from input)')
parser.add_option('-t', '--type', dest='type',
                  help='data type of array (default: float64)',
                  default='float64')
parser.add_option('-c', '--compress', dest='do_compress', action='store_true',
                  help='compress the dataset')
(opts, args) = parser.parse_args()

if len(args) != 3:
    sys.stderr.write('exactly 3 arguments are required\n')
    parser.print_help(sys.stderr)
    sys.exit(2)
    
try:
    dtype = getattr(numpy, opts.type)
except AttributeError:
    sys.stderr.write('invalid data type (%r) specified\n' % opts.type)
    sys.exit(2)
    
input_filename = args[0]
hdf5_filename = args[1]
hdf5_nodename = args[2]

if input_filename == '-':
    input_file = sys.stdin
else:
    input_file = open(input_filename, 'rt')

h5file = h5py.File(hdf5_filename)

ncol = opts.ncol
linebuffer = []
if ncol is None:
    line1 = input_file.readline()
    while line1 and line1[0] in ('#', '@'):
        line1 = input_file.readline()
        continue
    linebuffer.append(line1)
    line2 = input_file.readline()
    linebuffer.append(line2)
    fields1 = line1.split()
    fields2 = line2.split()
    if ncol is None:
        ncol = len(fields1)
        sys.stdout.write('using %d columns \n' % ncol)

if ncol < 10:
    chunk_size = 32768
else:
    chunk_size = 16384
conv_buffer = numpy.empty((chunk_size,ncol), dtype)
chunks_read = 0
lines_read = 0

nodenames = hdf5_nodename.split('/')
grp = None
for nodename in nodenames[:-1]:
    grp = (grp or h5file).require_group(nodename)

dset = (grp or h5file).create_dataset(nodenames[-1], (1,ncol), dtype,
                                      maxshape=(None, ncol),
                                      chunks=(chunk_size,ncol),
                                      compression = (opts.do_compress and 9 or None),
                                      )
for line in chain(linebuffer, input_file):
    if line[0] in ('#', '@'): continue
    
    fields = line.split()
    lines_read += 1
    if ncol == 1:
        conv_buffer[lines_read-1,0] = dtype(fields[0])
    else:
        conv_buffer[lines_read-1,:] = [dtype(field) for field in fields[0:ncol]]
    if lines_read == chunk_size:
        new_shape = (chunk_size*(chunks_read+1),) + tuple(dset.shape[1:])
        dset.resize(new_shape)
        dset[-chunk_size:,:] = conv_buffer
        lines_read = 0
        chunks_read += 1
        
if lines_read:
    old_shape = dset.shape
    new_shape = (dset.shape[0] + lines_read,) + dset.shape[1:]
    dset.resize(new_shape)
    dset[-lines_read:,:] = conv_buffer[:lines_read]
        
h5file.close()
input_file.close()
    

