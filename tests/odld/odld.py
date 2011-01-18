from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy, scipy
import wemd

class ODLDPropagator(wemd.propagators.WEMDPropagator):
    def __init__(self):
        self.variance = 0.001
        self.sigma = self.variance**0.5
        self.nsteps = 10
        
        self.ndim = 1
        
        self.potential_centers = None
        self.potential_widths = None
        self.potential_heights = None
        
    def propagate_segments(self, segments):
        if not segments:
            return
        assert len(segments[0].pcoord) == 1
        
        nsegs = len(segments)
        
        new_pcoords = numpy.empty((self.nsteps+1, nsegs, self.ndim), numpy.float64)
        
        for (iseg, segment) in enumerate(segments):
            new_pcoords[0,iseg,:] = segment.pcoord[-1]
        
        for istep in xrange(1, self.nsteps+1):
            random_displacements = numpy.random.normal(scale = self.sigma, size = (nsegs, self.ndim))
            new_pcoords[istep,:,:] = new_pcoords[istep-1,:,:] + random_displacements
        
        for (iseg, segment) in enumerate(segments):
            segment.pcoord = new_pcoords[:,iseg,:]

def construct_regions():
    from wemd.pcoords import ParticleSet, RectilinearRegionSet
    
    
    

def main():
    import itertools
    #sim_manager = wemd.sim_manager.WESimManager()
    
    # Construct driver tree
    
    # Construct simulation
    
    # Run simulation
    
    prop = ODLDPropagator()
    
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


if __name__ == '__main__':
    main()
