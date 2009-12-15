import sys, re
import numpy

class ArrayAccumulator(object):
    def __init__(self, dtype = None, buffer_size = 1024):
        self.data = None
        self.dtype = dtype
        self.buffer = list()
        self.buffer_size = max(0, buffer_size) or sys.maxint
        
    def append(self, fields):
        self.buffer.append(fields)
        if len(self.buffer) >= self.buffer_size:
            self.merge_buffer()
            
    def extend(self, iterable):
        self.buffer.extend(iterable)
        if len(self.buffer) >= self.buffer_size:
            self.merge_buffer()
            
    def merge_buffer(self):
        if not self.buffer: return
        try:
            oldshape = self.data.shape
        except AttributeError:
            self.data = numpy.array(self.buffer, self.dtype)
        else:
            self.data.resize((oldshape[0]+len(self.buffer), oldshape[1]))
            self.data[-len(self.buffer):,...] = self.buffer
        self.buffer = list()
        
class TextArrayConverter(object):
    re_valid_line_area = re.compile(r'^\s*([+-]?\d[^#]*)')
    re_field_sep = re.compile(r'\s+')
    
    def convert(self, stream, dtype = None, row_buffer_size = 1024):
        nlines = 0
        nfields = 0
        lineno = 0

        re_valid_line_area = self.re_valid_line_area
        re_field_sep = self.re_field_sep
        
        acc = ArrayAccumulator(dtype, row_buffer_size)
        for line in stream:
            lineno += 1
            m = re_valid_line_area.match(line)
            if m:
                nlines += 1
                fields = re_field_sep.split(m.group(1).strip())
                if nfields and nfields != len(fields):
                    raise TypeError('size mismatch at line %d' % lineno)
                else:
                    nfields = len(fields)
                acc.append(fields)
        acc.merge_buffer()
        return acc.data

if __name__ == '__main__':
    import os
    try:
        import cPickle as pickle
    except:
        import pickle
        
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [TXTFILE] [OUTFILE]',
                          description = 'convert a 2-D matrix stored as text'
                          +' into a binary format (pickled numpy array or '
                          +'raw binary data with a header).')
    parser.add_option('-b', '--buffer', dest='buffer_size', type='int',
                      help='read BUFFER_SIZE lines at a time (use 0 to read '
                           +'entire file before converting to an array)')
    parser.add_option('-t', '--dtype', dest='dtype',
                      help='use DTYPE as the data type of the numpy array '
                          +'(see numpy documentation for valid values)')
    
    (opts, args) = parser.parse_args()
    
    if opts.dtype:
        try:
            dtype = getattr(numpy, opts.dtype)
        except AttributeError:
            sys.stderr.write('invalid data type specified\n')
            sys.exit(2)
    else:
        dtype = numpy.float_
    
    if opts.buffer_size is not None:
        buffer_size = opts.buffer_size
    else:
        buffer_size = 1024    
    
    pickle_protocol = pickle.HIGHEST_PROTOCOL
    if len(args) > 2:
        parser.print_help(sys.stderr)
        sys.exit(2)
    elif len(args) in (1,2):
        infilename = args[0]
        try:
            outfilename = args[1]
        except IndexError:
            outfilename = '%s.dat' % (os.path.splitext(infilename)[0],)
        if infilename == '-':
            infile = sys.stdin
        else:
            infile = open(infilename, 'rt')
            
        if outfilename == '-':
            outfile = sys.stdout
        else:
            outfile = open(outfilename, 'wb')
    else: # no args
        infile = sys.stdin
        outfile = sys.stdout
        
    if outfile.isatty():
        sys.stderr.write('notice: writing pickle in text format because '
                         +'stdout is a tty\n')
        pickle_protocol = 0
        
    ac = TextArrayConverter()
    data = ac.convert(infile, dtype, buffer_size)
    
    pickle.dump(data, outfile, pickle_protocol)
    sys.exit(0)
