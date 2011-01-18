from __future__ import division; __metaclass__ = type

import numpy, operator

weight_getter = numpy.frompyfunc(operator.attrgetter('weight'), 1, 1)
count_getter = numpy.frompyfunc(operator.attrgetter('count'), 1, 1)

class ParticleSet(set):
    def __init__(self, iterable = None):
        if iterable is not None:
            super(ParticleSet,self).__init__(iterable)

    def map_to_bins(self, coords):
        # Degenerate case to terminate recursive descent
        return [self] * len(coords)

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
        return sum(region.count() for region in self.regions)        
    
    def transform_coords(self, coords):
        """Transform input coordinates into the coordinates spanned by this RegionSet. By default, does
        nothing (returns coords unchanged).  Note that to override this function, one may either subclass
        and override directly, or simply assign a function with signature transform(self, coords) to 
        transform_coords.  Also note that this could be as simple as returning a slice from coords 
        (i.e. to use only a subset of the total dimensionality of the coordinate)."""
        return coords
    
    def map_coords(self, coords):
        """ 
        Map each coord in coords to an index into self.regions identifying which region to which the corresponding
        coord belongs
        """
        raise NotImplementedError
        
    def map_to_bins(self, coords):
        """
         Recursively descend each region in this set until each coord in coords is mapped to a ParticleSet
         In other words, calling map_to_bins(coords) returns a list of the same length as coords where
         each entry is a ParticleSet; calling map_coords on a top-level RegionSet returns the bins to which 
         the particles belong.
        """
        
        bins = numpy.empty((len(coords),), numpy.object_)
        region_indices = self.map_coords(coords)
        
        for icoord in xrange(0, len(bins)):
            region = self.regions[region_indices[icoord]]
            # Descend recursively; since map_to_bins maps iterables to iterables, we must pass
            # a vector of length one and dereference the returned vector of length 1
            bins[icoord] = region.map_to_bins(coords[icoord:icoord+1])[0]
        
        if (bins == None).any():
            raise ValueError('one or more coords could not be mapped to bins')


class RectilinearRegionSet(RegionSet):
    """A multidimensional RegionSet divided by rectilinear boundaries"""
    
    def __init__(self):
        # A numpy array of RegionSets/ParticleSets
        self.region_array = None
        
        # A numpy array of integers
        self.target_counts_array = None
        
        self.ndim = None
        
        # A list of lists of bin boundaries
        # First dimension: dimension of coord
        # Remaining dimensions: bin boundaries
        # e.g. for 2-d boundaries at -1, 0, 1:
        # [ [-1, 0, 1], [-1, 0, 1] ]
        self.boundaries = None
        
        # Indirection array.  Each entry of this ndarray is an index into self.region_array.flat
        self.indir = None
                
    def construct_regions(self, boundaries):
        self.ndim = len(boundaries)
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
    
    def map_coords(self, coords):
        coords = numpy.asarray(self.transform_coords(coords))
        assert coords.ndim == 2
        assert coords.shape[1] == self.ndim
        region_indices = numpy.empty(coords.shape, numpy.uintp)
        flat_indices = numpy.empty(len(coords,), numpy.uintp)
                
        for idim in xrange(0, self.ndim):
            region_indices[:,idim] = numpy.digitize(coords[:,idim], self.boundaries[idim])
            if ( (region_indices[:, idim] == len(self.boundaries[idim]))
                |(region_indices[:, idim] == 0) ).any():
                # Beyond the upper bin limit
                raise ValueError('coord outside of bin space in dimension %d' % idim)
            region_indices[:,idim] -= 1
    
        # We now have an array of n-dimensional indices into our bin space, with each row corresponding to
        # one entry in coords
        indir = self.indir
        for icoord in xrange(0, len(coords)):
            flat_indices[icoord] = indir[region_indices[icoord]]
        
        return flat_indices
