from __future__ import division, print_function; __metaclass__ = type

import numpy, operator
import logging
log = logging.getLogger(__name__)

class ParticleCollection:
    def __init__(self, target_count = None, label = None):
        # How many particles should live here
        self.target_count = target_count
        
        # A descriptive label
        self.label = label
        
    def map_to_bins(self, coords):
        # Degenerate case to terminate recursive descent
        return [self] * len(coords)
    
    def get_all_bins(self):
        return [self]

    def reweight(self, new_weight):
        """Reweight all particles in this collection so that the total weight is new_weight"""
        current_weight = self.weight
        for p in self:
            p.weight *= new_weight / current_weight
    
    @property
    def weight(self):
        'Total weight of all particles in this collection'
        weight = 0.0
        for particle in self:
            weight += particle.weight
        return weight
    
    @weight.setter
    def weight(self, new_weight):
        return self.reweight(new_weight)
        
    @property
    def count(self):
        """The number of particles in this collection"""
        return len(self)
    
    def __repr__(self):
        return '<{classname} at 0x{id:x}, label={label!r}>'.format(classname=self.__class__.__name__, 
                                                                   id=id(self), 
                                                                   label=self.label)
    
class ParticleSet(set, ParticleCollection):
    def __init__(self, iterable = None, target_count = None, label = None):
        ParticleCollection.__init__(self, target_count, label)
        if iterable is not None:
            set.__init__(self, iterable)
        
class ParticleList(list, ParticleCollection):
    def __init__(self, iterable = None, target_count = None, label = None):
        ParticleCollection.__init__(self, target_count, label)
        if iterable is not None:
            list.__init__(self, iterable)
        
    def clear(self):
        self[:] = []
                
class RegionSet:    
    def __init__(self):
        
        # Regions (RegionSets or ParticleSets)
        self.regions = None
        
        # Transform coordinates to those used in this RegionSet; could be as simple as
        # returning a slice, for systems where nested coordinate regions handle different
        # parts of the pcoord data as stored in the HDF5 file.
        self.transform_coords = None
        
        # Used to select fast variants of mapping algorithms when topologies
        # aren't nested
        self._simple_topology = False
        self._stale_topology = True
        self._binmap = None
        
    def clear(self):
        for region in self.regions:
            region.clear()
                                        
    def reweight(self, new_weight):
        """Reweight all particles in this set so that the total weight is new_weight"""
        current_weight = self.weight
        for r in self.regions:
            r.weight *= new_weight / current_weight

    @property
    def weight(self):
        'Total weight of all particles in this set' 
        return sum(region.weight for region in self.regions)
    
    @weight.setter
    def weight(self, new_weight):
        return self.reweight(new_weight)
        
    @property
    def count(self):
        """The number of particles in this set"""
        return sum(region.count for region in self.regions)        
        
    def prep_coords(self, coords):
        """Prepare a list of coordinates for further processing. Ensures that the coordinates
        are filtered through self.transform_coords(), and further that input and output 
        conforms to assumptions about array shape -- namely, that an array of coords is
        two-dimensional, with the first dimension indexing time and the second dimension
        indexing progress coordinate dimension."""

        # ensure we have an array        
        coords = numpy.asarray(coords)
        
        # ensure we have a 2-d array
        if coords.ndim == 1:
            coords = numpy.expand_dims(coords, axis=1)
        elif coords.ndim > 2:
            raise TypeError('coords must be a 1- or 2-d array')
        
        # apply transform_coords if it's defined
        try:
            transform_coords = self.transform_coords
        except AttributeError:
            # No transformation defined
            return coords
        else:
            if transform_coords is None:
                # Another way of expressing "no transformation defined"
                return coords
        
        # a transformation is defined; apply it
        coords = numpy.asarray(transform_coords(coords))
        
        # and again ensure that we have a 2-d array to return
        if coords.ndim == 1:
            coords = numpy.expand_dims(coords, axis=1)
        elif coords.ndim > 2:
            raise TypeError('transform_coords() must return a 1- or 2-d array')
        return coords
    
    def scan_topology(self):
        """Scan the topology of this bin space. Must be called prior to calling the various map_to_...
        functions"""
        
        all_bins = self.get_all_bins()
        if list(all_bins) == list(self.regions):
            self._simple_topology = True
            self._binmap = None
            log.info('using algorithms optimized for non-nested bin spaces')
        else:
            self._simple_topology = False
            self._binmap = {id(bin): ibin for (ibin,bin) in enumerate(all_bins)}
            log.info('using algorithms appropriate for nested bin spaces')
        self._stale_topology = False
        
    def _map_to_indices(self, coords):
        """ 
        Map each coord in coords to an index into self.regions identifying the region to which each corresponding
        coord belongs.  The first index denotes separate sets of coordinates, and the second index denotes 
        a component within a coordinate set.  Subclasses (must) override this to implement binning functionality. 
        """
        raise NotImplementedError
          
    def map_to_bins(self, coords):
        """
        Recursively descend each region in this set until each coord in coords is mapped to a ParticleSet
        In other words, calling map_to_bins(coords) returns a list of the same length as coords where
        each entry is a ParticleSet; calling map_to_bins on a top-level RegionSet returns the bins to which 
        the particles belong.
        """
        
        if self._stale_topology:
            self.scan_topology()
        
        coords = self.prep_coords(coords)
        bins = numpy.empty((len(coords),), numpy.object_)
        region_indices = self._map_to_indices(coords)
        
        if self._simple_topology:
            bins[:] = [self.regions[index] for index in region_indices]
        else:
            populated_indices = set(region_indices)
            for index in populated_indices:
                # create a vector of len(bins) of booleans indicating which rows we're updating 
                bin_selector = (region_indices == index)
                bins[bin_selector] = self.regions[index].map_to_bins(coords[bin_selector])
        
        if None in bins:
            raise ValueError('one or more coords could not be mapped to bins')
        
        return bins

    def get_all_bins(self):
        """Get a list of all bins contained in this RegionSet, descending depth-first"""
        bins = []
        for region in self.regions:
            bins.extend(region.get_all_bins())
        return bins        

    def map_to_all_indices(self, coords):
        """Map sets of coordinates to indices into the list of bins returned by `get_all_bins()`
        
        This function caches a map from bin identity to index. If the topology of regions changes
        between calls to map_to_all_indices(), scan_topology() must also be called to ensure that
        the map gets updated properly.
        
        Currently, there are no optimizations for non-nested bin topologies. Such optimizations may
        increase throughput for analysis of relatively simple topologies."""

        if self._stale_topology:
            self.scan_topology()
            
        if self._simple_topology:
            # No nesting, so we can delegate to _map_to_indices()
            return self._map_to_indices(self.prep_coords(coords))
        else:
            return [self._binmap[id(bin)] for bin in self.map_to_bins(coords)]
                
    def get_bin_containing(self, coord):
        return self.map_to_bins(self.prep_coords([coord]))[0]
    
    def replace_region(self, coord, new_container):
        """Deprecated synonym for replace_region_containing()"""
        import warnings
        warnings.warn('replace_region() is deprecated in favor of replace_region_containing()', DeprecationWarning,2)
        self.replace_region_containing(coord, new_container)
        
    def replace_region_containing(self, coord, new_container):
        """Replace the region containing coord with new_container"""
        idx = self._map_to_indices(self.prep_coords([coord]))[0]
        self.regions[idx] = new_container
        self._stale_topology = True

class PiecewiseRegionSet(RegionSet):
    """A multidimensional RegionSet which uses a set of functions to map coordinates to regions.
    In the event multiple functions match for a given coordinate set, the first one matched (which in turn
    is the first one passed to the constructor) takes precedence."""
    
    def __init__(self, functions, regions = None, container_class = ParticleSet):
        super(PiecewiseRegionSet,self).__init__()
        self.functions = list(functions)
        self.regions = numpy.empty((len(functions),), numpy.object_)
        if regions is not None:
            self.regions[:] = regions
        for (ireg, region) in enumerate(self.regions):
            if region is None:
                self.regions[ireg] = container_class()            

    def _map_to_indices(self, coords):
        assert coords.ndim == 2
        functions = self.functions
        
        region_indices = numpy.empty((len(coords),), numpy.uintp)
        
        for (icoord, coord) in enumerate(coords):
            for (ifunc, func) in enumerate(functions):
                if func(coord):
                    region_indices[icoord] = ifunc
                    break
            else:
                # loop terminated without break, which means no function matched
                raise ValueError('coordinate ({!r}) outside of bin space'.format(coord))

        return region_indices
        
class RectilinearRegionSet(RegionSet):
    """A multidimensional RegionSet divided by rectilinear (interior) boundaries"""
    
    def __init__(self, boundaries, container_class = ParticleSet):
        # A numpy array of RegionSets/ParticleSets
        self.region_array = None
        
        # because of the accessor for region, region_array must be set before
        # calling our parent's __init__
        super(RectilinearRegionSet,self).__init__()
                   
        self.ndim = None
        
        # A list of lists of bin boundaries
        # First dimension: dimension of coord
        # Remaining dimensions: bin boundaries
        # e.g. for 2-d boundaries at -1, 0, 1:
        # [ [-1, 0, 1], [-1, 0, 1] ]
        self.boundaries = None
                
        if boundaries is not None:
            self.construct_regions(boundaries, container_class)
                                
    def construct_regions(self, boundaries, container_class):
        self.ndim = len(boundaries)
        self.region_array = numpy.empty(tuple(len(boundary_entry)-1 for boundary_entry in boundaries), numpy.object_)
        self.boundaries = [numpy.array(boundary_set) for boundary_set in boundaries]
        
        for index in numpy.ndindex(self.region_array.shape):
            bin = container_class()

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
            
        # flat iterators, like self.regions, always index in C order
        # (create a numpy Fortran-order array and call flat on it and watch
        # what happens -- only if you iterate on flat in C order will the
        # data come out right)
        # use clever cumulative products based on shape dotted with the
        # indices to reduce N-d indices to flat indices         
        self._extents = numpy.ones((self.region_array.ndim,), numpy.uintp)
        self._extents[:-1] = numpy.cumprod(self.region_array.shape[::-1])[:-1][::-1]
                    
    # The mandated region and target_counts iterables are 1-dimensional, so provide
    # accessors to one-dimensional representations
    @property
    def regions(self):
        return self.region_array.flat
    
    @regions.setter
    def regions(self, regions):
        if regions is not None:
            self.region_array.flat = regions
    
    def _map_to_indices(self, coords):
        assert coords.ndim == 2

        if coords.shape[1] != self.ndim:
            raise TypeError('length of coordinate tuples ({}) does not match the dimensionality of this region set ({})'
                            .format(coords.shape[1], self.ndim))

        flat_indices   = numpy.zeros((len(coords),), numpy.uintp)
        dim_indices = numpy.empty((len(coords),), numpy.uintp)
        extents = self._extents
        
        for idim in xrange(0, self.ndim):
            dim_indices[:] = numpy.digitize(coords[:,idim], self.boundaries[idim])
            if ( (dim_indices == len(self.boundaries[idim])) | (dim_indices == 0) ).any():
                # Beyond the bin limits
                raise ValueError('coordinate outside of bin space in dimension %d' % idim)
            # Adjust for how numpy.digitize chooses to return its values
            dim_indices -= 1
            flat_indices += dim_indices*extents[idim]
                        
        return flat_indices
