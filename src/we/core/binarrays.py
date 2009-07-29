import itertools
import numpy
from particles import ParticleCollection
__metaclass__ = type

class Bin(ParticleCollection):
    def __init__(self, iterable = None, 
                 bin_id = None, index = None,
                 ideal_num = None, split_threshold = None,
                 merge_threshold_min = None, merge_threshold_max = None):
        super(Bin,self).__init__(iterable or [])
        self.bin_id = bin_id
        self.index = index
        self.ideal_num = ideal_num
        self.split_threshold = split_threshold
        self.merge_threshold_min = merge_threshold_min
        self.merge_threshold_max = merge_threshold_max
        
    def __repr__(self):
        return '<%s(%s) index=%r, %d particles, norm=%g>' \
               % (self.__class__.__name__, hex(id(self)), self.index,
                  len(self), self.norm)

class BinArray:
    def __init__(self, boundaries, ideal_num, 
                 split_threshold, merge_threshold_min, merge_threshold_max):
        self.boundaries = boundaries
        self.shape = tuple(bound.shape[-1]-1 for bound in boundaries)
        self.ndim = len(self.shape)
        
        self.bins = bins = numpy.empty([bound.shape[-1]-1 for bound in boundaries],
                                       dtype = numpy.object_)
        
        self.ideal_num = self._per_bin(ideal_num)
        self.split_threshold = self._per_bin(split_threshold)
        self.merge_threshold_min = self._per_bin(merge_threshold_min)
        self.merge_threshold_max = self._per_bin(merge_threshold_max)
        
        for i in xrange(0, len(bins.flat)):
            index = numpy.unravel_index(i, bins.shape)
            bins[index] = Bin([], bin_id = i, index=index,
                              ideal_num = self.ideal_num[index],
                              split_threshold = self.split_threshold[index],
                              merge_threshold_min = self.merge_threshold_min[index],
                              merge_threshold_max = self.merge_threshold_max[index])
            
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
        
    def population_array(self):
        pop_array = numpy.empty(self.bins.shape, numpy.float_)
        for (i, bin) in enumerate(self.bins.flat):
            pop_array.flat[i] = bin.norm
        return pop_array
            
    def nparticles_array(self):
        nparticles_array = numpy.empty(self.bins.shape, numpy.int_)
        for (i, bin) in enumerate(self.bins.flat):
            nparticles_array.flat[i] = len(bin)
            
    def clear(self):
        for bin in self.bins.flat:
            bin.clear()

    def __iter__(self):
        return self.bins.flat