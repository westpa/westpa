from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy, scipy
import wemd

class ODLDPropagator(wemd.propagators.WEMDPropagator):
    def __init__(self):
        self.variance = 0.001
        self.sigma = self.variance**0.5
        self.nsteps = 10
        
        self.ndim = 0
        
        self.potential_centers = None
        self.potential_widths = None
        self.potential_heights = None
        
    def propagate_segments(self, segments):
        if not segments:
            return
        assert len(segments[0].pcoord) == 1
        
        nsegs = len(segments)
        
        new_pcoords = numpy.empty((self.nsteps+1, nsegs, self.ndim), numpy.float64)
        new_pcoords[0,:,:] = [segment.pcoord for segment in segments]
        
        for istep in xrange(1, self.nsteps+1):
            random_displacements = numpy.random.normal(scale = self.sigma, size = (nsegs, self.ndim))
            new_pcoords[istep,:,:] = new_pcoords[istep-1,:,:] + random_displacements
        
        for (iseg, segment) in enumerate(segments):
            segment.pcoord = new_pcoords[:,iseg,:] = new_pcoords[:,iseg,:]

def construct_regions():
    from wemd.pcoords import ParticleSet, RectilinearRegionSet
    
    
    

#def main():
    #sim_manager = wemd.sim_manager.WESimManager()
    
    # Construct driver tree
    
    # Construct simulation
    
    # Run simulation
    