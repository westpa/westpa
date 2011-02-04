from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy, scipy
import wemd
from wemd import Segment, WEMDSystem
from wemd.pcoords import PiecewiseRegionSet, ParticleSet, RectilinearRegionSet

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

class ODLDPropagator(wemd.propagators.WEMDPropagator):
    def __init__(self, runtime_config, sim_manager = None):
        self.runtime_config = runtime_config
        self.sim_manager = sim_manager
        
        self.variance = 0.001
        self.sigma = self.variance**0.5
        self.nsteps = 10
        
        self.ndim = 1
        
        self.potential_centers = None
        self.potential_widths = None
        self.potential_heights = None
        
    def propagate(self, segments):
        if not segments:
            return
        
        nsegs = len(segments)
        
        new_pcoords = numpy.empty((self.nsteps+1, nsegs, self.ndim), numpy.float64)
        
        for (iseg, segment) in enumerate(segments):
            new_pcoords[0,iseg,:] = segment.pcoord[0]
        
        for istep in xrange(1, self.nsteps+1):
            random_displacements = numpy.random.normal(scale = self.sigma, size = (nsegs, self.ndim))
            new_pcoords[istep,:,:] = new_pcoords[istep-1,:,:] + random_displacements
        
        for (iseg, segment) in enumerate(segments):
            segment.pcoord = new_pcoords[:,iseg,:]
            segment.status = Segment.SEG_STATUS_COMPLETE

class ODLDSystem(WEMDSystem):
    def __init__(self, sim_manager):
        self.sim_manager = sim_manager
        
        # Construct bin regions
        # Divide space into x<0 and x>=0
        self.region_set = PiecewiseRegionSet([(lambda x: x<0), (lambda x: x>=0)])
        
        # Replace each bin with a set of smaller bins
        self.region_set.replace_region([-0.5], RectilinearRegionSet([[float('-inf')] + list(numpy.arange(-1.0, 0.1, 0.1))]))
        self.region_set.replace_region([0.5], RectilinearRegionSet(boundaries=[list(numpy.arange(0.0, 1.1, 0.1)) + [float('inf')]]))
        
        # 10 particles per bin
        for bin in self.region_set.get_all_bins():
            bin.target_count = 20
            
        # Except for target bins, which are denoted by a target_count of zero; here, abs(x) > 1
        left_target = self.region_set.get_bin_containing([-1.1])
        right_target = self.region_set.get_bin_containing([1.1])
        left_target.target_count = right_target.target_count = 0
        left_target.label = 'left_target'
        right_target.label = 'right_target'
        
        self.target_states = [('left_target', left_target),
                              ('right_target', right_target)]
        
        self.initial_states = [('left', 0.9, [-0.1], self.region_set.get_bin_containing([-0.1])),
                                     ('right', 0.1, [0.1], self.region_set.get_bin_containing([0.1]))]
        self.pcoord_ndim = 1
        self.pcoord_len = 11
        self.pcoord_dtype = numpy.float64
        


def test_regions_simple(segments):
    rrs = RectilinearRegionSet([numpy.arange(-1.0,1.1,0.1)])
    
    print(repr(rrs.boundaries))
    print(repr(rrs.region_array))
    
    for segment in segments:
        pci = rrs.map_to_indices([segment.pcoord[-1,:]])
        print('segment', segment.seg_id, 'pcoord', segment.pcoord[-1], 'bin', pci)
        
    bins = rrs.map_to_bins([segment.pcoord[-1] for segment in segments])
    for (bin,segment) in zip(bins,segments):
        bin.add(segment)
    
    print(repr(rrs.region_array))

def test_regions_nested(segments): 
    rs = PiecewiseRegionSet([(lambda x: x<0), (lambda x: x>=0)])
    rs.replace_region([-0.5], RectilinearRegionSet(boundaries=[numpy.arange(-1.0, 0.1, 0.1)]))
    rs.replace_region([0.5], RectilinearRegionSet(boundaries=[numpy.arange(0.0, 1.1, 0.1)]))
                
    bins = rs.map_to_bins([segment.pcoord[-1] for segment in segments])
    for (bin,segment) in zip(bins,segments):
        bin.add(segment)
    
    print(repr(rs.regions))
    for (ireg, region) in enumerate(rs.regions):
        print('region', ireg, 'count', region.count)
        for (iireg, iregion) in enumerate(region.regions):
            print('iregion', iireg, 'count', iregion.count)


def test_nopot():
    import itertools
    
    prop = ODLDPropagator({})
    
    parent = None
    seg_id_generator = itertools.count(1)
    n_segs = 4
    segs = [wemd.Segment(seg_id=seg_id_generator.next(),p_parent=parent,pcoord=numpy.zeros((1,1), numpy.float64))
                for iseg in xrange(0, n_segs)]

    for n in xrange(0, 10):
        segs = [wemd.Segment(seg_id=seg_id_generator.next(),p_parent=parent,pcoord=seg.pcoord[-1])
                    for seg in segs]        
        prop.propagate_segments(segs)
        print('n =', n, 'pcoords:')
        numpy.savetxt(sys.stdout, numpy.column_stack([seg.pcoord for seg in segs]))
        parent = seg
        pcoord = seg.pcoord[-1]
    
    return segs

def run_basic_tests():
    segs = test_nopot()
    test_regions_simple(segs)
    test_regions_nested(segs)

if __name__ == '__main__':
    run_basic_tests()