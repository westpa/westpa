from __future__ import division; __metaclass__ = type

import numpy, operator

weight_getter = numpy.frompyfunc(operator.attrgetter('weight'), 1, 1)
count_getter = numpy.frompyfunc(operator.attrgetter('count'), 1, 1)

class ParticleSet(set):
    def __init__(self, iterable = None):
        if iterable is not None:
            super(ParticleSet,self).__init__(iterable)

    def map_to_bins(self, pcoords):
        # Degenerate case to terminate recursive descent
        return [self] * len(pcoords)

    def reweight(self, new_weight):
        """Reweight all particles in this set so that the total weight is new_weight"""
        current_weight = self.weight
        for p in self:
            p.weight *= new_weight / current_weight
    
    @property
    def weight(self):
        'Total weight of all particles in this set'
        return numpy.add.reduce(weight_getter(self))
    
    weight.setter = reweight
        
    @property
    def count(self):
        """The number of particles in this set"""
        return len(self)

            
class RegionSet:    
    def __init__(self):
        # Regions (RegionSets or ParticleSets)
        self.regions = None
        
        # Target counts for each of the above regions
        # An entry of None denotes no target count to be enforced
        self.target_counts = None

    def reweight(self, new_weight):
        """Reweight all particles in this set so that the total weight is new_weight"""
        current_weight = self.weight
        for r in self.regions:
            r.weight *= new_weight / current_weight
    
    @property
    def weight(self):
        'Total weight of all particles in this set' 
        return numpy.add.reduce(map(weight_getter, self.regions))
    
    weight.setter = reweight
        
    @property
    def count(self):
        """The number of particles in this set"""
        return len(self)        
    
    def construct_regions(self, boundaries):
        raise NotImplementedError
    
    def map_pcoords(self, pcoords):
        """ 
        Map each pcoord in pcoords to an index into self.regions identifying which region to which the corresponding
        pcoord belongs
        """
        raise NotImplementedError
        
    def map_to_bins(self, pcoords):
        """
         Recursively descend each region in this set until each pcoord in pcoords is mapped to a ParticleSet
         In other words, calling map_to_bins(pcoords) returns a list of the same length as pcoords where
         each entry is a ParticleSet; calling map_pcoords on a top-level RegionSet returns the bins to which 
         the particles belong.
        """
        
        bins = numpy.empty((len(pcoords),), numpy.object_)
        region_indices = self.map_pcoords(pcoords)
        
        for ipcoord in xrange(0, len(bins)):
            region = self.regions[region_indices[ipcoord]]
            # Descend recursively; since map_to_bins maps iterables to iterables, we must pass
            # a vector of length one and dereference the returned vector of length 1
            bins[ipcoord] = region.map_to_bins(pcoords[ipcoord:ipcoord+1])[0]
        
        if (bins == None).any():
            raise ValueError('one or more pcoords could not be mapped to bins')


class RectilinearRegionSet(RegionSet):
    """A multidimensional RegionSet divided by rectilinear boundaries"""
    
    def __init__(self):
        # A numpy array of RegionSets/ParticleSets
        self.region_array = None
        
        # A numpy array of integers
        self.target_counts_array = None
        
        # A list of lists of bin boundaries
        # First dimension: dimension of pcoord
        # Remaining dimensions: bin boundaries
        # e.g. for 2-d boundaries at -1, 0, 1:
        # [ [-1, 0, 1], [-1, 0, 1] ]
        self.boundaries = None
        
        # Indirection array.  Each entry of this ndarray is an index into self.region_array.flat
        self.indir = None
        
    def construct_regions(self, boundaries):
        self.region_array = numpy.empty(tuple(len(boundary_entry)-1 for boundary_entry in boundaries), numpy.object_)
        self.target_counts_array = numpy.zeros(self.region_array.shape, numpy.uint)
        
        for index in numpy.ndindex(self.region_array.shape):
            self.region_array[index] = ParticleSet()
            
        # Admittedly not the most memory-efficient, but at least it's quite simple
        self.indir = numpy.arange(self.region_array.size).reshape(self.region_array.shape)
                            
    # The mandated region and target_counts iterables are 1-dimensional, so provide
    # accessors to one-dimensional representations
    @property
    def regions(self):
        return self.region_array.flat
    
    @property
    def target_counts(self):
        return self.target_counts_array.flat
    
    def map_pcoords(self, pcoords):
        pcoords = numpy.asarray(pcoords)
        assert pcoords.ndim == 2
        region_indices = numpy.empty(pcoords.shape, numpy.uintp)
        flat_indices = numpy.empty(len(pcoords,), numpy.uintp)
                
        for idim in xrange(0, self.boundaries.ndim):
            region_indices[:,idim] = numpy.digitize(pcoords[:,idim], self.boundaries[idim])
            if ( (region_indices[:, idim] == len(self.boundaries[idim]))
                |(region_indices[:, idim] == 0) ).any():
                # Beyond the upper bin limit
                raise ValueError('pcoord outside of bin space in dimension %d' % idim)
            region_indices[:,idim] -= 1
    
        # We now have an array of n-dimensional indices into our bin space, with each row corresponding to
        # one entry in pcoords
        indir = self.indir
        for ipcoord in xrange(0, len(pcoords)):
            flat_indices[ipcoord] = indir[region_indices[ipcoord]]
        
        return flat_indices
    
                    
            