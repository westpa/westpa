import itertools
import numpy
from particles import ParticleCollection
__metaclass__ = type

class Bin(ParticleCollection):    
    def __init__(self, iterable = None, 
                 bin_id = None, index = None,
                 ideal_num = None):
        super(Bin,self).__init__(iterable or [])
        self.bin_id = bin_id
        self.index = index
        self.ideal_num = ideal_num
        
    def __repr__(self):
        return '<%s(%s) index=%r, %d particles, norm=%g>' \
               % (self.__class__.__name__, hex(id(self)), 
                  self.index,
                  len(self), self.norm)

class BinArray:
    def __init__(self, boundaries, ideal_num):
        self.boundaries = boundaries
        self.shape = tuple(bound.shape[-1]-1 for bound in boundaries)
        self.ndim = len(self.shape)
        
        self.bins = bins = numpy.empty([bound.shape[-1]-1 for bound in boundaries],
                                       dtype = numpy.object_)
        
        self.ideal_num = self._per_bin(ideal_num)
        
        for i in xrange(0, len(bins.flat)):
            index = numpy.unravel_index(i, bins.shape)
            bins[index] = Bin([], bin_id = i, index=index,
                              ideal_num = self.ideal_num[index])
            
    def _per_bin(self, const_or_ndarray):
        o = const_or_ndarray
        try:
            if o.shape == self.bins.shape:
                return o
            else:
                raise TypeError('array shapes do not match')
        except AttributeError:
            nd = numpy.empty(self.bins.shape, type(o))
            nd.fill(o)
            return nd
        
    def get_norm(self):
        return sum(b.norm for b in self.bins)
    
    norm = property(get_norm, None, None, 'total weight in all bins')
    
    def renorm(self, weight=1.0):
        wfactor = weight/self.get_norm()
        for bin in self:
            for particle in bin:
                particle.weight *= wfactor
        
    def population_array(self):
        pop_array = numpy.empty(self.bins.shape, numpy.float_)
        for (i, bin) in enumerate(self.bins.flat):
            pop_array.flat[i] = bin.norm
        return pop_array
            
    def nparticles_array(self):
        nparticles_array = numpy.empty(self.bins.shape, numpy.int_)
        for (i, bin) in enumerate(self.bins.flat):
            nparticles_array.flat[i] = len(bin)
        return nparticles_array
    
    def clear(self):
        for bin in self.bins.flat:
            bin.clear()
            
    def __len__(self):
        return self.bins.size

    def __iter__(self):
        return self.bins.flat