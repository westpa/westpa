import xdrlib
import numpy

class GenericArrayData(object):
    MAGIC = 'wemdGAD.0001'
    HEADER_SIZE = 40
    
    def __init__(self):
        self.dtype = None
        self.ndim = None
        self.data_offset = None
        self.shape = None
        self.data = None
        
    def calc_data_offset(self, shape):
        ndim = len(shape)
        return self.HEADER_SIZE + 8*ndim
        
    def write_header_to(self, stream, t0 = None, dt = None, 
                        dtype = None, shape = None):
        ndim = self.ndim = len(shape)
            
        dtype = dtype or self.dtype
        if dtype is None:
            dtypestr = ''
        else:
            dtypestr = numpy.dtype(dtype).descr[0][1]

        packer = xdrlib.Packer()
        packer.pack_fstring(16, self.MAGIC)
        packer.pack_fstring(8, dtypestr)
        packer.pack_uhyper(ndim)
        self.data_offset = self.calc_data_offset(shape)
        packer.pack_uhyper(self.data_offset)        
        packer.pack_farray(ndim, shape, packer.pack_uhyper)
        
        packed = packer.get_buffer()
        assert len(packed) == self.data_offset
        
        stream.write(packed)
        
    def read_header_from(self, stream):
        unpacker = xdrlib.Unpacker(stream.read(self.HEADER_SIZE))
        magic = unpacker.unpack_fstring(16)
        if magic[0:len(self.MAGIC)] != self.MAGIC:
            raise ValueError('invalid magic number %r' % magic)
        dtypestr = unpacker.unpack_fstring(8)
        dtypestr = dtypestr.strip('\x00')
        if not dtypestr:
            self.dtype = None
        else:
            self.dtype = numpy.dtype(dtypestr)
        self.ndim = unpacker.unpack_uhyper()
        self.data_offset = unpacker.unpack_uhyper()
        unpacker.done()
        unpacker = xdrlib.Unpacker(stream.read(self.ndim * 8))
        self.shape = tuple(unpacker.unpack_farray(self.ndim, unpacker.unpack_uhyper))
        unpacker.done()
        
    def mmap_array_from(self, stream, mode='r'):
        if self.dtype is None:
            self.read_header_from(stream)
        
        return numpy.memmap(stream.name, mode='r',
                            dtype = self.dtype, shape = self.shape,
                            offset = self.data_offset)
            
class UniformTimeData(object):
    MAGIC = 'wemdUTD.0001'
    HEADER_SIZE = 56
    
    def __init__(self):
        self.dtype = None
        self.ndim = None
        self.data_offset = None
        self.t0 = None
        self.dt = None
        self.shape = None
        self.data = None
        
    def calc_data_offset(self, shape):
        ndim = len(shape)
        return self.HEADER_SIZE + 8*ndim
        
    def write_header_to(self, stream, dtype = None, shape = None):
        ndim = self.ndim = len(shape)
            
        dtype = dtype or self.dtype
        if dtype is None:
            dtypestr = ''
        else:
            dtypestr = numpy.dtype(dtype).descr[0][1]

        packer = xdrlib.Packer()
        packer.pack_fstring(16, self.MAGIC)
        packer.pack_fstring(8, dtypestr)
        packer.pack_uhyper(ndim)
        self.data_offset = self.calc_data_offset(shape)
        packer.pack_uhyper(self.data_offset)        
        packer.pack_farray(ndim, shape, packer.pack_uhyper)
        
        packed = packer.get_buffer()
        assert len(packed) == self.data_offset
        
        stream.write(packed)
        
    def read_header_from(self, stream):
        unpacker = xdrlib.Unpacker(stream.read(self.HEADER_SIZE))
        magic = unpacker.unpack_fstring(16)
        if magic[0:len(self.MAGIC)] != self.MAGIC:
            raise ValueError('invalid magic number %r' % magic)
        dtypestr = unpacker.unpack_fstring(8)
        dtypestr = dtypestr.strip('\x00')
        if not dtypestr:
            self.dtype = None
        else:
            self.dtype = numpy.dtype(dtypestr)
        self.t0 = unpacker.unpack_double()
        self.dt = unpacker.unpack_double()
        self.ndim = unpacker.unpack_uhyper()
        self.data_offset = unpacker.unpack_uhyper()
        unpacker.done()
        unpacker = xdrlib.Unpacker(stream.read(self.ndim * 8))
        self.shape = tuple(unpacker.unpack_farray(self.ndim, unpacker.unpack_uhyper))
        unpacker.done()
        
    def mmap_array_from(self, stream, mode='r'):
        if self.dtype is None:
            self.read_header_from(stream)
        
        return numpy.memmap(stream.name, mode='r',
                            dtype = self.dtype, shape = self.shape,
                            offset = self.data_offset)
        
    
        
        
        
        
        
    
        
        
        
        
        
        
        