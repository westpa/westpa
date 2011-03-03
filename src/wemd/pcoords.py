from __future__ import division, print_function; __metaclass__ = type

import numpy, operator
import logging
log = logging.getLogger(__name__)

weight_getter = numpy.frompyfunc(operator.attrgetter('weight'), 1, 1)
count_getter = numpy.frompyfunc(operator.attrgetter('count'), 1, 1)

class ParticleSet(set):
    def __init__(self, iterable = None, target_count = None, label = None):
        if iterable is not None:
            super(ParticleSet,self).__init__(iterable)
        
        # How many particles should live here
        self.target_count = target_count
        
        # A descriptive label
        self.label = label
        
    def map_to_bins(self, coords):
        # Degenerate case to terminate recursive descent
        coords = list(coords)
        return [self] * len(coords)
    
    def get_all_bins(self):
        return [self]

    def reweight(self, new_weight):
        """Reweight all particles in this set so that the total weight is new_weight"""
        current_weight = self.weight
        for p in self:
            p.weight *= new_weight / current_weight
    
    @property
    def weight(self):
        'Total weight of all particles in this set'
        #return numpy.add.reduce(map(weight_getter, self))
        weight = 0.0
        for particle in self:
            weight += particle.weight
        return weight
    
    @weight.setter
    def weight(self, new_weight):
        return self.reweight(new_weight)
        
    @property
    def count(self):
        """The number of particles in this set"""
        return len(self)
    
    def __repr__(self):
        return '<ParticleSet at 0x{id:x}, label={label!r}>'.format(id=id(self), label=self.label)
    
            
class RegionSet:    
    def __init__(self):
        
        # Regions (RegionSets or ParticleSets)
        self.regions = None
    
    def clear(self):
        for region in self.regions:
            region.clear()
        
    def get_all_bins(self):
        bins = []
        for region in self.regions:
            bins.extend(region.get_all_bins())
        return bins
        
                                
    def reweight(self, new_weight):
        """Reweight all particles in this set so that the total weight is new_weight"""
        current_weight = self.weight
        for r in self.regions:
            r.weight *= new_weight / current_weight

    @property
    def weight(self):
        'Total weight of all particles in this set' 
        return numpy.add.reduce(map(weight_getter, self.regions))
    
    @weight.setter
    def weight(self, new_weight):
        return self.reweight(new_weight)
        
    @property
    def count(self):
        """The number of particles in this set"""
        return sum(region.count for region in self.regions)        
    
    def transform_coords(self, coords):
        """Transform input coordinates into the coordinates spanned by this RegionSet. By default, does
        nothing (returns coords unchanged).  Note that to override this function, one may either subclass
        and override directly, or simply assign a function with signature transform(self, coords) to 
        transform_coords.  Also note that this could be as simple as returning a slice from coords 
        (i.e. to use only a subset of the total dimensionality of the coordinate)."""
        return coords
    
    def map_to_indices(self, coords):
        """ 
        Map each coord in coords to an index into self.regions identifying which region to which the corresponding
        coord belongs.  The first index denotes separate sets of coordinates, and the second index denotes 
        a component within a coordinate set.
        """
        raise NotImplementedError
        
    def map_to_bins(self, coords):
        """
         Recursively descend each region in this set until each coord in coords is mapped to a ParticleSet
         In other words, calling map_to_bins(coords) returns a list of the same length as coords where
         each entry is a ParticleSet; calling map_to_indices on a top-level RegionSet returns the bins to which 
         the particles belong.
        """
        #coords = numpy.array(coords)
        coords = list(coords)
        bins = numpy.empty((len(coords),), numpy.object_)
        region_indices = self.map_to_indices(coords)
        
        for icoord in xrange(0, len(bins)):
            region = self.regions[region_indices[icoord]]
            # Descend recursively; since map_to_bins maps iterables to iterables, we must pass
            # a vector of length one and dereference the returned vector of length 1
            bins[icoord] = region.map_to_bins(coords[icoord:icoord+1])[0]
        
        if None in bins:
            raise ValueError('one or more coords could not be mapped to bins')
        
        return bins
    
    def get_bin_containing(self, coord):
        return self.map_to_bins([coord])[0]
    
    def replace_region(self, coord, new_container):
        """Replace the region containing coord with new_container"""
        idx = self.map_to_indices([coord])[0]
        self.regions[idx] = new_container

class PiecewiseRegionSet(RegionSet):
    """A multidimensional RegionSet which uses a set of functions to map coordinates to regions"""
    def __init__(self, functions, regions = None):
        super(PiecewiseRegionSet,self).__init__()
        self.functions = list(functions)
        self.regions = numpy.empty((len(functions),), numpy.object_)
        if regions is not None:
            self.regions[:] = regions
        for (ireg, region) in enumerate(self.regions):
            if region is None:
                self.regions[ireg] = ParticleSet()            
        
    def map_to_indices(self, coords):
        coords = numpy.asarray(self.transform_coords(coords))
        region_pred = numpy.empty((len(coords), len(self.functions)), numpy.bool_)
        
        for (ifunc, func) in enumerate(self.functions):
            rsl = numpy.apply_along_axis(func, 1, coords)
            region_pred[:,ifunc] = rsl.flat
        
        # why does this seem too fancy?
        apred = numpy.argwhere(region_pred)
        if not (apred[:,0] == numpy.arange(0, len(coords))).all():
            raise ValueError('coordinate outside of bin space')
        region_indices = apred[:,1]
        return region_indices
        
class RectilinearRegionSet(RegionSet):
    """A multidimensional RegionSet divided by rectilinear (interior) boundaries"""
    
    def __init__(self, boundaries):
        
        # A numpy array of RegionSets/ParticleSets
        self.region_array = None
                
        self.ndim = None
        
        # A list of lists of bin boundaries
        # First dimension: dimension of coord
        # Remaining dimensions: bin boundaries
        # e.g. for 2-d boundaries at -1, 0, 1:
        # [ [-1, 0, 1], [-1, 0, 1] ]
        self.boundaries = None
        
        # Indirection array.  Each entry of this ndarray is an index into self.region_array.flat
        self.indir = None
        
        if boundaries is not None:
            self.construct_regions(boundaries)
                                
    def construct_regions(self, boundaries):
        self.ndim = len(boundaries)
        self.region_array = numpy.empty(tuple(len(boundary_entry)-1 for boundary_entry in boundaries), numpy.object_)
        self.boundaries = [numpy.array(boundary_set) for boundary_set in boundaries]
        
        for index in numpy.ndindex(self.region_array.shape):
            bin = ParticleSet()

            # This list comprehension reduces to this:
            # index is (i_dim0, i_dim1, idim_2, i_dim3, ... )
            #bounds=[]
            #for idim in xrange(0,len(index)):
            #    lb = boundaries[idim][index[idim]]
            #    ub = boundaries[idim][index[idim]+1]
            #    bounds.append[(lb,ub)]
            bounds = [(boundaries[idim][index[idim]], boundaries[idim][index[idim]+1])
                      for idim in xrange(0,self.ndim)]
            bin.label = repr(bounds)
            self.region_array[index] = bin
            
        # Admittedly not the most memory-efficient, but at least it's quite simple and fast
        self.indir = numpy.arange(self.region_array.size,dtype=numpy.uintp).reshape(self.region_array.shape)
        
    # The mandated region and target_counts iterables are 1-dimensional, so provide
    # accessors to one-dimensional representations
    @property
    def regions(self):
        return self.region_array.flat
    
    @regions.setter
    def regions(self, regions):
        self.region_array.flat = regions
    
    def map_to_indices(self, coords):
        coords = numpy.atleast_2d(self.transform_coords(coords))
        assert coords.ndim == 2
        assert coords.shape[1] == self.ndim
        region_indices = numpy.empty(coords.shape, numpy.uintp)
        flat_indices   = numpy.zeros((len(coords),), numpy.uintp)
        
        # flat iterators, like self.regions, always index in C order
        # this sets up the ordering
        extents = numpy.ones((self.region_array.ndim,), numpy.uintp)
        extents[:-1] = numpy.cumprod(self.region_array.shape[::-1])[:-1][::-1]
        
        for idim in xrange(0, self.ndim):
            region_indices[:,idim] = numpy.digitize(coords[:,idim], self.boundaries[idim])
            if ( (region_indices[:, idim] == len(self.boundaries[idim]))
                |(region_indices[:, idim] == 0) ).any():
                # Beyond the bin limits
                raise ValueError('coordinate outside of bin space in dimension %d' % idim)
            
            # Adjust for how numpy.digitize chooses to return its values
            region_indices[:,idim] -= 1
            flat_indices += region_indices[:,idim]*extents[idim]
                        
        return flat_indices
        
if __name__ == '__main__':
    regions = RectilinearRegionSet([[-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2]])
    
    coords = -4*numpy.random.random(size=(10,3)) + 2
    indices = regions.map_to_indices(coords)
    
    for (coord, index) in zip(coords,indices):
        print('coords:', coord, 'index:', index, 'bin:', regions.regions[index])
    